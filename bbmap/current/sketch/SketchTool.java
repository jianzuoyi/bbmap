package sketch;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.atomic.AtomicInteger;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.ByteStreamWriter;
import kmer.HashArray1D;
import kmer.HashForest;
import kmer.KmerNode;
import kmer.KmerTableSet;
import shared.Parser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ByteBuilder;
import stream.KillSwitch;
import structures.LongHeap;
import structures.LongList;

/**
 * @author Brian Bushnell
 * @date June 28, 2016
 *
 */
public final class SketchTool extends SketchObject {
	
	/*--------------------------------------------------------------*/
	/*----------------         Main Method          ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		
		args=Parser.parseConfig(args);
		if(Parser.parseHelp(args, true)){
			//printOptions();
			System.exit(0);
		}
		
		Timer t=new Timer();
		t.start();
		
		//Create a new CountKmersExact instance
		SketchTool mhs=new SketchTool(args);
		t.stop();
		System.err.println("Time: \t"+t);
	}
	
	//For testing
	public SketchTool(String[] args){
		System.err.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");
		
		/* Set global defaults */
		ReadWrite.ZIPLEVEL=2;
		ReadWrite.USE_UNPIGZ=true;
		
		/* Initialize local variables with defaults */
		Parser parser=new Parser();
		
		ArrayList<String> list=new ArrayList<String>();
		
		int size_=SketchObject.targetSketchSize;
		int minKeyOccuranceCount_=1;
		float cutoff=0.02f;
		int mode=ONE_SKETCH;
		
		/* Parse arguments */
		for(int i=0; i<args.length; i++){

			final String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if("null".equalsIgnoreCase(b)){b=null;}
			while(a.charAt(0)=='-' && (a.indexOf('.')<0 || i>1 || !new File(a).exists())){a=a.substring(1);}
			
			if(a.equals("in")){
				if(b!=null){
					for(String s : b.split(",")){
						list.add(s);
					}
				}
			}else if(a.equals("minkmercount") || a.equals("minkmeycount")){
				minKeyOccuranceCount_=Integer.parseInt(b);
			}else if(a.equals("cutoff")){
				cutoff=Float.parseFloat(b);
			}else if(parseSketchFlags(arg, a, b)){
				//Do nothing
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else if(b==null){
				list.add(arg);
			}else{
				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		stTargetSketchSize=size_;
		minKeyOccuranceCount=minKeyOccuranceCount_;
		
		Timer t=new Timer();
		ArrayList<Sketch> sketches=loadSketches_MT(mode, 1, -1, list);
		t.stop();
		System.err.println("Loaded "+sketches.size()+" sketches in \t"+t);
		t.start();
		Sketch sketch=sketches.get(0);
		CompareBuffer buffer=new CompareBuffer(false);
		for(int i=1; i<sketches.size(); i++){
			Sketch sketch2=sketches.get(i);
			final int matches=sketch.countMatches(sketch2, buffer, null, false);
			assert(matches==buffer.matches);
			
			if(matches>=3){
				final float idWeighted=matches/(float)buffer.minDivisor();
				final float idMin=matches/(float)buffer.minSize();
				final float idMax=matches/(float)buffer.maxSize();
				
				if(idWeighted>=cutoff && matches>=3){
					System.out.println(String.format(Locale.ROOT, "%.2f%% WID, %.2f%% MINID, %.2f%% MAXID, %d matches, %d length",
							100*idWeighted, 100*idMin, 100*idMax, matches, buffer.minDivisor())+" for "+sketch.name()+" vs "+sketch2.name());
				}
			}
		}
		t.stop();
		System.err.println("Compared "+(sketches.size()-1)+" sketches in \t"+t);
		
//		String fname=list.get(0);
//		Sketch sketch=loadSketches(fname).get(0);
//		for(int i=1; i<list.size(); i++){
//			String fname2=list.get(i);
//			Sketch sketch2=loadSketches(fname2).get(0);
//			float identity=sketch.identity(sketch2);
//			System.out.println("Identity for "+fname+" vs "+fname2+":\t"+String.format(Locale.ROOT, "%.2f%%", 100*identity));
//		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------      Normal Constructor      ----------------*/
	/*--------------------------------------------------------------*/
	
	public SketchTool(int size_, int minKeyOccuranceCount_){
		stTargetSketchSize=size_;
		minKeyOccuranceCount=minKeyOccuranceCount_;
		
		assert(!amino || !rcomp) : "rcomp should be false in amino mode.";
		assert(!amino || (k*AminoAcid.AMINO_SHIFT<64)) : "Protein sketches require 1 <= K <= "+(63/AminoAcid.AMINO_SHIFT)+".";
		assert(k>0 && k<=32) : "Sketches require 1 <= K <= 32."; //123
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/
	
	public Sketch toSketch(KmerTableSet tables, boolean multithreaded){
		final int threads=(multithreaded ? Tools.mid(1, Shared.threads(), tables.ways()) : 1);
		return (threads<2 ? toSketch_ST(tables) : toSketch_MT(tables, threads));
	}
	
	private Sketch toSketch_ST(KmerTableSet tables){
		SketchHeap heap=(stTargetSketchSize>0 ? new SketchHeap(stTargetSketchSize, minKeyOccuranceCount) : null);
		LongList list=new LongList();

		KmerTableSet kts=(KmerTableSet)tables;
		for(int tnum=0; tnum<kts.ways; tnum++){
			HashArray1D table=kts.getTable(tnum);
			if(stTargetSketchSize>0){
				toHeap(table, heap);
			}else{
				toList(table, list);
			}
		}
		return stTargetSketchSize>0 ? new Sketch(heap, false) : toSketch(list);
	}
	
	private Sketch toSketch_MT(KmerTableSet tables, final int threads){
		ArrayList<SketchThread> alst=new ArrayList<SketchThread>(threads);
		AtomicInteger ai=new AtomicInteger(0);
		for(int i=0; i<threads; i++){
			alst.add(new SketchThread(ai, tables));
		}

		//Start the threads
		for(SketchThread pt : alst){
			pt.start();
		}

		ArrayList<SketchHeap> heaps=new ArrayList<SketchHeap>(threads);
		LongList list=new LongList();
		
		for(SketchThread pt : alst){

			//Wait until this thread has terminated
			while(pt.getState()!=Thread.State.TERMINATED){
				try {
					//Attempt a join operation
					pt.join();
				} catch (InterruptedException e) {
					//Potentially handle this, if it is expected to occur
					e.printStackTrace();
				}
			}
			
			if(stTargetSketchSize>=0){
				if(pt.heap!=null && pt.heap.size()>0){
					heaps.add(pt.heap);
				}
			}else{
				if(pt.list!=null){list.append(pt.list);}
				pt.list=null;
			}
		}
		alst.clear();
		
		return stTargetSketchSize>=0 ? toSketch(heaps) : toSketch(list);
	}
	
	public static final long[] toSketchArray(SketchHeap sheap){
		int maxLen=toSketchSize(sheap.genomeSizeBases, sheap.genomeSizeKmers, sheap.genomeSizeEstimate(), SketchObject.targetSketchSize);
//		System.err.println(sheap.genomeSizeBases+", "+sheap.genomeSizeKmers+", "+sheap.genomeSizeEstimate()+", "+SketchObject.size+" -> "+maxLen); //123
		return toSketchArray(sheap, maxLen);
	}
	
	public static final long[] toSketchArrayOld(SketchHeap sheap, int maxLen){
		final LongHeap heap=sheap.heap;
		final int initial=heap.size();
		final int len=Tools.min(maxLen, initial);
		final long[] array=KillSwitch.allocLong1D(len);
		
		int toSkip=heap.size()-len;
		for(int i=0; i<toSkip; i++){heap.poll();}
		for(int i=0; i<len; i++){
			array[i]=Long.MAX_VALUE-heap.poll();
		}
		Tools.reverseInPlace(array);
		assert(heap.size()==0) : heap.size()+", "+len+", "+maxLen+", "+initial;
		return array;
	}
	
	public static final long[] toSketchArray(SketchHeap sheap, int maxLen){
		if(sheap.setMode){return toSketchArrayOld(sheap, maxLen);}
		long[] keys=sheap.map().toArray(sheap.minKeyOccuranceCount);
		for(int i=0; i<keys.length; i++){
//			assert(keys[i]>0) : Arrays.toString(keys);
			keys[i]=Long.MAX_VALUE-keys[i];
//			assert(keys[i]>0) : Arrays.toString(keys);
		}
		Arrays.sort(keys);
		if(keys.length>maxLen){
			keys=Arrays.copyOf(keys, maxLen);
		}
		
		final LongHeap heap=sheap.heap;
		heap.clear();
		assert(heap.size()==0) : heap.size()+", "+maxLen;
		return keys;
	}
	
	public SketchHeap toHeap(HashArray1D table, SketchHeap heap){
//		if(heap==null){heap=new LongHeap(size, true);}
		long[] kmers=table.array();
		int[] counts=table.values();
		for(int i=0; i<table.arrayLength(); i++){
			int count=counts[i];
			if(count>=minKeyOccuranceCount){
				heap.genomeSizeKmers++;
				long hash=hash(kmers[i]);
				if(hash>=minHashValue){
					heap.add(hash);
				}
			}
		}
		HashForest forest=table.victims();
		if(forest!=null){
			for(KmerNode kn : forest.array()){
				if(kn!=null){addRecursive(heap, kn);}
			}
		}
		return heap;
	}
	
	public LongList toList(HashArray1D table, LongList list){
//		if(heap==null){heap=new LongHeap(size, true);}
		long[] kmers=table.array();
		int[] counts=table.values();
		for(int i=0; i<table.arrayLength(); i++){
			int count=counts[i];
			if(count>=minKeyOccuranceCount){
				long hash=hash(kmers[i]);
				if(hash>=minHashValue){
					list.add(hash);
				}
			}
		}
		HashForest forest=table.victims();
		if(forest!=null){
			for(KmerNode kn : forest.array()){
				if(kn!=null){addRecursive(list, kn);}
			}
		}
		return list;
	}
	
//	public long[] toSketchArray(ArrayList<LongHeap> heaps){
//		if(heaps.size()==1){return toSketchArray(heaps.get(0));}
//		LongList list=new LongList(size);
//		for(LongHeap heap : heaps){
//			while(heap.size()>0){list.add(Long.MAX_VALUE-heap.poll());}
//		}
//		list.sort();
//		list.shrinkToUnique();
//		list.size=Tools.min(size, list.size);
//		return list.toArray();
//	}
	
	public Sketch toSketch(ArrayList<SketchHeap> heaps){
		SketchHeap a=heaps.get(0);
		for(int i=1; i<heaps.size(); i++){
			SketchHeap b=heaps.get(i);
			a.add(b);
		}
//		assert(false) : a.size()+", "+new Sketch(a).array.length;
		return new Sketch(a, false);
	}
	
	public Sketch toSketch(LongList list){
		list.sort();
		assert(list.size==0 || list.get(list.size()-1)>=minHashValue) : list.size+", "+list.get(list.size()-1)+", "+minHashValue;
		list.shrinkToUnique();
		list.reverse();
		for(int i=0; i<list.size; i++){list.array[i]=Long.MAX_VALUE-list.array[i];}
		return new Sketch(list.toArray());
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Helpers            ----------------*/
	/*--------------------------------------------------------------*/
	
	private void addRecursive(SketchHeap heap, KmerNode kn){
		if(kn==null){return;}
		if(kn.count()>=minKeyOccuranceCount){
			heap.genomeSizeKmers++;
			long hash=hash(kn.pivot());
			if(hash>=minHashValue){heap.add(hash);}
		}
		if(kn.left()!=null){addRecursive(heap, kn.left());}
		if(kn.right()!=null){addRecursive(heap, kn.right());}
	}
	
	private void addRecursive(LongList list, KmerNode kn){
		if(kn==null){return;}
		if(kn.count()>=minKeyOccuranceCount){
			long hash=hash(kn.pivot());
			if(hash>=minHashValue){list.add(hash);}
		}
		if(kn.left()!=null){addRecursive(list, kn.left());}
		if(kn.right()!=null){addRecursive(list, kn.right());}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------             I/O              ----------------*/
	/*--------------------------------------------------------------*/
	
//	public static ArrayList<Sketch> loadSketches_ST(String...fnames){
//		ArrayList<Sketch> sketches=null;
//		for(String s : fnames){
//			ArrayList<Sketch> temp;
//			if(s.indexOf(',')<0 || s.startsWith("stdin") || new File(s).exists()){
//				temp=loadSketches(s);
//			}else{
//				temp=loadSketches_ST(s.split(","));
//			}
//			if(sketches==null){sketches=temp;}
//			else{sketches.addAll(temp);}
//		}
//		return sketches;
//	}
	
//	public static ArrayList<Sketch> loadSketches_MT(ArrayList<String> fnames){
//		return loadSketches_MT(0, null, fnames.toArray(new String[0]));
//	}
	
	public ArrayList<Sketch> loadSketches_MT(int mode, float samplerate, long reads, ArrayList<String> fnames){
		return loadSketches_MT(mode, samplerate, reads, fnames.toArray(new String[0]));
	}
	
	public ArrayList<Sketch> loadSketches_MT(int mode, float samplerate, long reads, String...fnames){
		ConcurrentLinkedQueue<String> decomposedFnames=new ConcurrentLinkedQueue<String>();
		for(String s : fnames){
			if(s.indexOf(',')<0 || s.startsWith("stdin") || new File(s).exists()){
				decomposedFnames.add(s);
			}else{
				for(String s2 : s.split(",")){
					decomposedFnames.add(s2);
				}
			}
		}

		if(decomposedFnames.size()==0){return null;}
		if(decomposedFnames.size()==1){return loadSketches(decomposedFnames.poll(), null, mode, samplerate, reads);}
		
		
		//Determine how many threads may be used
		final int threads=Tools.min(Shared.threads(), decomposedFnames.size());

		//Fill a list with LoadThreads
		ArrayList<LoadThread> allt=new ArrayList<LoadThread>(threads);
		
		for(int i=0; i<threads; i++){
			allt.add(new LoadThread(decomposedFnames, mode, samplerate, reads));
		}
		
		ArrayList<Sketch> sketches=new ArrayList<Sketch>();
		
		//Start the threads
		for(LoadThread lt : allt){lt.start();}

		//Wait for completion of all threads
		boolean success=true;
		for(LoadThread lt : allt){

			//Wait until this thread has terminated
			while(lt.getState()!=Thread.State.TERMINATED){
				try {
					//Attempt a join operation
					lt.join();
				} catch (InterruptedException e) {
					//Potentially handle this, if it is expected to occur
					e.printStackTrace();
				}
			}
			sketches.addAll(lt.list);
			success&=lt.success;
		}
		assert(success) : "Failure loading some files.";
		return sketches;
	}
	
	private class LoadThread extends Thread{
		
		public LoadThread(ConcurrentLinkedQueue<String> queue_, int mode_, float samplerate_, long reads_) {
			queue=queue_;
			list=new ArrayList<Sketch>();
//			map=map_;
//			if(map==null){list=new ArrayList<Sketch>();}
			smm=new SketchMakerMini(SketchTool.this, mode_);
			samplerate=samplerate_;
			reads=reads_;
		}
		
		public void run(){
			success=false;
			for(String fname=queue.poll(); fname!=null; fname=queue.poll()){
				ArrayList<Sketch> temp=null;
				try {
					temp=loadSketches(fname, smm, smm.mode, samplerate, reads);
				} catch (Throwable e) {
					System.err.println("Failure loading "+fname+":\n"+e);
					e.printStackTrace();
					success=false;
				}
				if(temp!=null){
					for(Sketch s : temp){add(s);}
				}
			}
			success=true;
		}
		
		private void add(Sketch s){
			if(list!=null){
				list.add(s);
				return;
			}
			assert(false) : "Unsupported."; //The map logic is broken; needs to be synchronized.
//			if(s.taxID<0){return;}
////			assert(s.taxID>-1) : s.toHeader();
//			TaxNode tn=tree.getNode(s.taxID);
//			while(tn!=null && tn.pid!=tn.id && tn.level<taxLevel){
//				TaxNode temp=tree.getNode(tn.pid);
//				if(temp==null){break;}
//				tn=temp;
//			}
//			if(tn==null){return;}
//			Integer key=tn.id;
//			Sketch old=map.get(key);
//			if(old==null){
//				s.taxID=key;
//				map.put(key, s);
//			}else{
//				synchronized(old){
//					old.add(s, maxLen);
//				}
//			}
		}
		
		final ConcurrentLinkedQueue<String> queue;
		ArrayList<Sketch> list;
		boolean success=false;
		final SketchMakerMini smm;
		final float samplerate;
		final long reads;
		
//		ConcurrentHashMap<Integer, Sketch> map;
		
	}
	
	public ArrayList<Sketch> loadSketches(String fname, SketchMakerMini smm, int mode, float samplerate, long reads){
		FileFormat ff=FileFormat.testInput(fname, FileFormat.FASTA, null, false, true);
		
		if(ff.fasta() || ff.fastq() || ff.samOrBam()){
			if(smm==null){smm=new SketchMakerMini(this, mode);}
			ArrayList<Sketch> sketches=smm.toSketches(fname, samplerate, reads);
			return sketches;
		}
		
		boolean A48=Sketch.CODING==Sketch.A48, HEX=Sketch.CODING==Sketch.HEX, delta=Sketch.delta;
		
		ArrayList<Sketch> sketches=new ArrayList<Sketch>();
		ByteFile bf=ByteFile.makeByteFile(fname, false);
		int currentSketchSize=stTargetSketchSize;
		int taxID=-1;
		long spid=-1;
		long imgID=-1;
		long genomeSizeBases=0, genomeSizeKmers=0, genomeSequences=0;
		String name=null, name0=null;
		LongList list=null;
		long sum=0;
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length>0){
//				System.err.println("Processing line "+new String(line));
				if(line[0]=='#'){
					if(list!=null){
						assert(list.size==list.array.length);
						list.shrink();
						if(list.size()>0){
							Sketch sketch=new Sketch(list.array, taxID, imgID, genomeSizeBases, genomeSizeKmers, genomeSequences, name, name0, ff.simpleName());
							sketches.add(sketch);
						}
//						System.err.println("Made sketch "+sketch);
					}
					name=name0=null;
					list=null;
					sum=0;
					taxID=-1;
					imgID=-1;
					genomeSizeBases=0;
					genomeSizeKmers=0;
					int k_sketch=defaultK;
					int k2_sketch=0;
					
					if(line.length>1){
						String[] split=new String(line, 1, line.length-1).split("\t");
						for(String s : split){
							final String sub=s.substring(s.indexOf(':')+1);
							if(s.startsWith("SZ:") || s.startsWith("SIZE:")){//Sketch length
								currentSketchSize=Integer.parseInt(sub);
							}else if(s.startsWith("SZ:")){
								currentSketchSize=Integer.parseInt(sub);
							}else if(s.startsWith("CD:")){//Coding
								A48=HEX=delta=false;
								
								for(int i=0; i<sub.length(); i++){
									char c=sub.charAt(i);
									if(c=='A'){A48=true;}
									else if(c=='H'){HEX=true;}
									else if(c=='R'){A48=HEX=false;}
									else if(c=='D'){delta=true;}
									else if(c=='M'){assert(amino) : "Amino sketch in non-amino mode: "+new String(line);}
									else if(c=='8'){assert(amino8) : "Amino8 sketch in non-amino8 mode: "+new String(line);}
									else{assert(false) : "Unknown coding symbol: "+c+"\t"+new String(line);}
								}
								
							}else if(s.startsWith("K:")){//Kmer length
								if(sub.indexOf(',')>=0){
									String[] subsplit=sub.split(",");
									assert(subsplit.length==2) : "Bad header component "+s;
									int x=Integer.parseInt(subsplit[0]);
									int y=Integer.parseInt(subsplit[1]);
									k_sketch=Tools.max(x, y);
									k2_sketch=Tools.min(x, y);
								}else{
									k_sketch=Integer.parseInt(sub);
									k2_sketch=0;
								}
							}else if(s.startsWith("H:")){//Hash version
								int hashVersion_sketch=Integer.parseInt(sub);
								assert(hashVersion_sketch==HASH_VERSION) : 
									"Sketch hash_version "+hashVersion_sketch+" differs from loaded hash_version "+HASH_VERSION+"\n"+new String(line);
							}else if(s.startsWith("GS:") || s.startsWith("GSIZE:")){//Genomic bases
								genomeSizeBases=Long.parseLong(sub);
							}else if(s.startsWith("GK:") || s.startsWith("GKMERS:")){//Genomic kmers
								genomeSizeKmers=Long.parseLong(sub);
							}else if(s.startsWith("GQ:")){
								genomeSequences=Long.parseLong(sub);
							}else if(s.startsWith("GE:")){//Genome size estimate kmers
								//ignore
							}else if(s.startsWith("ID:") || s.startsWith("TAXID:")){
								taxID=Integer.parseInt(sub);
							}else if(s.startsWith("IMG:")){
								imgID=Long.parseLong(sub);
							}else if(s.startsWith("SPID:")){
								spid=Integer.parseInt(sub);
							}else if(s.startsWith("NM:") || s.startsWith("NAME:")){
								name=sub;
							}else if(s.startsWith("FN:")){
								fname=sub;
							}else if(s.startsWith("NM0:")){
								name0=sub;
							}else{
								assert(false) : "Unsupported header tag "+s;
							}
						}
					}
					if(KILL_OK){
						if(k_sketch!=k){KillSwitch.kill("Sketch kmer length "+k_sketch+" differs from loaded kmer length "+k+"\n"+new String(line));}
						if(k2_sketch!=k2){KillSwitch.kill("Sketch kmer length "+k_sketch+","+k2_sketch+" differs from loaded kmer length "+k+","+k2+"\n"+new String(line));}
					}else{//Potential hang
						assert(k_sketch==k) : "Sketch kmer length "+k_sketch+" differs from loaded kmer length "+k+"\n"+new String(line);
						assert(k2_sketch==k2) : "Sketch kmer length "+k_sketch+","+k2_sketch+" differs from loaded kmer length "+k+","+k2+"\n"+new String(line);
					}
					if(currentSketchSize>0){list=new LongList(currentSketchSize);}
				}else{
					long x=(A48 ? Sketch.parseA48(line) : HEX ? Sketch.parseHex(line) : Tools.parseLong(line));
//					System.err.println("sum="+sum+", x="+x+" -> "+(sum+x));
					sum+=x;
					assert(x>=0) : x+"\n"+new String(line);
					assert(sum>=0 || !delta) : "The sketch was made with delta compression off.  Please regenerate it.";
					assert(list!=null) : new String(line);
					list.add(delta ? sum : x);
					//						System.err.println("List="+list);
				}
			}
		}
		if(list!=null && list.size>0){
			assert(list.size==list.array.length);
			list.shrink();
			Sketch sketch=new Sketch(list.array, taxID, imgID, genomeSizeBases, genomeSizeKmers, genomeSequences, name, name0, ff.simpleName());
			sketch.spid=spid;
			sketches.add(sketch);
		}
		return sketches;
	}
	
	public ArrayList<Sketch> loadSketchesFromString(String sketchString){
		boolean A48=Sketch.CODING==Sketch.A48, HEX=Sketch.CODING==Sketch.HEX, delta=Sketch.delta;
		
		ArrayList<Sketch> sketches=new ArrayList<Sketch>();
		int currentSketchSize=stTargetSketchSize;
		int taxID=-1;
		long spid=-1;
		long imgID=-1;
		long genomeSizeBases=0, genomeSizeKmers=0, genomeSequences=0;
		String name=null, name0=null, fname=null;
		LongList list=null;
		long sum=0;
		String[] split0=sketchString.split("\n");
		for(String line : split0){
			if(line.length()>0){
//				System.err.println("Processing line "+new String(line));
				if(line.charAt(0)=='#'){
					if(line.length()>1 && line.charAt(1)=='#'){
						//ignore
					}else{
						if(list!=null){
							assert(list.size==list.array.length);
							list.shrink();
							if(list.size()>0){
								Sketch sketch=new Sketch(list.array, taxID, imgID, genomeSizeBases, genomeSizeKmers, genomeSequences, name, name0, fname);
								sketches.add(sketch);
							}
							//						System.err.println("Made sketch "+sketch);
						}
						name=name0=null;
						list=null;
						sum=0;
						taxID=-1;
						imgID=-1;
						genomeSizeBases=0;
						genomeSizeKmers=0;
						genomeSequences=0;
						int k_sketch=defaultK;
						int k2_sketch=0;

						if(line.length()>1){
							String[] split=line.substring(1).split("\t");
							for(String s : split){
								final String sub=s.substring(s.indexOf(':')+1);
								if(s.startsWith("SZ:") || s.startsWith("SIZE:")){//Sketch length
									currentSketchSize=Integer.parseInt(sub);
								}else if(s.startsWith("SZ:")){
									currentSketchSize=Integer.parseInt(sub);
								}else if(s.startsWith("CD:")){//Coding
									A48=HEX=delta=false;
									
									for(int i=0; i<sub.length(); i++){
										char c=sub.charAt(i);
										if(c=='A'){A48=true;}
										else if(c=='H'){HEX=true;}
										else if(c=='R'){A48=HEX=false;}
										else if(c=='D'){delta=true;}
										else if(c=='M'){assert(amino) : "Amino sketch in non-amino mode: "+new String(line);}
										else if(c=='8'){assert(amino8) : "Amino8 sketch in non-amino8 mode: "+new String(line);}
										else{assert(false) : "Unknown coding symbol: "+c+"\t"+new String(line);}
									}
									
								}else if(s.startsWith("K:")){//Kmer length
									if(sub.indexOf(',')>=0){
										String[] subsplit=sub.split(",");
										assert(subsplit.length==2) : "Bad header component "+s;
										int x=Integer.parseInt(subsplit[0]);
										int y=Integer.parseInt(subsplit[1]);
										k_sketch=Tools.max(x, y);
										k2_sketch=Tools.min(x, y);
									}else{
										k_sketch=Integer.parseInt(s);
										k2_sketch=0;
									}
								}else if(s.startsWith("H:")){//Hash version
									int hashVersion_sketch=Integer.parseInt(sub);
									assert(hashVersion_sketch==HASH_VERSION) : 
										"Sketch hash_version "+hashVersion_sketch+" differs from loaded hash_version "+HASH_VERSION+"\n"+new String(line);
								}else if(s.startsWith("GS:") || s.startsWith("GSIZE:")){//Genomic bases
									genomeSizeBases=Long.parseLong(sub);
								}else if(s.startsWith("GK:") || s.startsWith("GKMERS:")){//Genomic kmers
									genomeSizeKmers=Long.parseLong(sub);
								}else if(s.startsWith("GQ:")){
									genomeSequences=Long.parseLong(sub);
								}else if(s.startsWith("GE:")){//Genome size estimate kmers
									//ignore
								}else if(s.startsWith("ID:") || s.startsWith("TAXID:")){
									taxID=Integer.parseInt(sub);
								}else if(s.startsWith("IMG:")){
									imgID=Long.parseLong(sub);
								}else if(s.startsWith("SPID:")){
									spid=Integer.parseInt(sub);
								}else if(s.startsWith("NM:") || s.startsWith("NAME:")){
									name=sub;
								}else if(s.startsWith("FN:")){
									fname=sub;
								}else if(s.startsWith("NM0:")){
									name0=sub;
								}else{
									assert(false) : "Unsupported header tag "+s;
								}
							}
						}
						if(KILL_OK){
							if(k_sketch!=k){KillSwitch.kill("Sketch kmer length "+k_sketch+" differs from loaded kmer length "+k+"\n"+new String(line));}
							if(k2_sketch!=k2){KillSwitch.kill("Sketch kmer length "+k_sketch+","+k2_sketch+" differs from loaded kmer length "+k+","+k2+"\n"+new String(line));}
						}else{//Potential hang
							assert(k_sketch==k) : "Sketch kmer length "+k_sketch+" differs from loaded kmer length "+k+"\n"+new String(line);
							assert(k2_sketch==k2) : "Sketch kmer length "+k_sketch+","+k2_sketch+" differs from loaded kmer length "+k+","+k2+"\n"+new String(line);
						}
						if(currentSketchSize>0){list=new LongList(currentSketchSize);}
					}
				}else{
					long x=(A48 ? Sketch.parseA48(line) : HEX ? Sketch.parseHex(line) : Long.parseLong(line));
//					System.err.println("sum="+sum+", x="+x+" -> "+(sum+x));
					sum+=x;
					assert(x>=0) : x+"\n"+new String(line);
					assert(sum>=0 || !delta) : "The sketch was made with delta compression off.  Please regenerate it.";
					list.add(delta ? sum : x);
					//						System.err.println("List="+list);
				}
			}
		}
		
		if(list!=null && list.size>0){
			assert(list.size==list.array.length || list.size()==0);
			list.shrink();
			Sketch sketch=new Sketch(list.array, taxID, imgID, genomeSizeBases, genomeSizeKmers, genomeSequences, name, name0, fname);
			sketch.spid=spid;
			sketches.add(sketch);
		}
		return sketches;
	}
	
	public static boolean write(ArrayList<Sketch> sketches, FileFormat ff[]){
		final int len=ff.length;
		ByteStreamWriter tsw[]=new ByteStreamWriter[len];
		for(int i=0; i<len; i++){
			tsw[i]=new ByteStreamWriter(ff[i]);
			tsw[i].start();
		}
		boolean error=false;
		for(int i=0; i<sketches.size(); i++){
			write(sketches.get(i), tsw[i%len]);
		}
		for(int i=0; i<len; i++){
			error|=tsw[i].poisonAndWait();
		}
		return error;
	}
	
	public static boolean write(ArrayList<Sketch> sketches, FileFormat ff){
		ByteStreamWriter tsw=new ByteStreamWriter(ff);
		tsw.start();
		for(Sketch sketch : sketches){
			write(sketch, tsw);
		}
		return tsw.poisonAndWait();
	}
	
	public static boolean write(Sketch sketch, FileFormat ff){
//		System.err.println(ff.name()+", "+new File(ff.name()).exists());
		ByteStreamWriter tsw=new ByteStreamWriter(ff);
//		assert(false) : new File(ff.name()).exists();
		tsw.start();
		write(sketch, tsw);
		return tsw.poisonAndWait();
	}
	
	public static void write(Sketch sketch, ByteStreamWriter tsw){
		write(sketch, tsw, new byte[12]);
	}
	
	public static void write(Sketch sketch, ByteStreamWriter tsw, byte[] temp){
		long prev=0;
		long[] array=sketch.array;
		final boolean A48=Sketch.CODING==Sketch.A48, HEX=Sketch.CODING==Sketch.HEX, delta=Sketch.delta;
		final ByteBuilder bb=sketch.toHeader().append('\n');
		for(int i=0; i<array.length; i++){
			long key=array[i];
			assert(key!=prev) : "Bad sketch - key "+key+" == "+prev;
			long x=key-prev;
			if(A48){
				Sketch.appendA48(x, bb, temp);
				bb.append('\n');
			}else if(HEX){
				bb.append(Long.toHexString(x));
				bb.append('\n');
			}else{
				bb.append(x);
				bb.append('\n');
			}
			if(delta){prev=key;}
		}
		tsw.print(bb);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Nested Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** Converts KmerTableSets to Heaps */
	private class SketchThread extends Thread {

		SketchThread(AtomicInteger next_, KmerTableSet kts_){
			next=next_;
			kts=kts_;
		}

		public void run(){
			final int ways=kts.ways();
			int tnum=next.getAndIncrement();
			while(tnum<ways){
				HashArray1D table=kts.getTable(tnum);
				if(stTargetSketchSize>0){
					if(heap==null){heap=new SketchHeap(stTargetSketchSize, minKeyOccuranceCount);}
					toHeap(table, heap);
				}else{
					if(list==null){list=new LongList();}
					toList(table, list);
				}
				tnum=next.getAndIncrement();
			}
		}

		final AtomicInteger next;
		final KmerTableSet kts;
		SketchHeap heap;
		LongList list;
	}
		
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	final int stTargetSketchSize;
	public final int minKeyOccuranceCount;
	
}
