package sketch;

import shared.Tools;
import stream.ByteBuilder;
import structures.AbstractBitSet;
import structures.LongList;
import tax.ImgRecord;

/**
 * @author Brian Bushnell
 * @date July 7, 2016
 *
 */
public class Sketch extends SketchObject implements Comparable<Sketch> {
	
	/*--------------------------------------------------------------*/
	/*----------------         Constructors         ----------------*/
	/*--------------------------------------------------------------*/
	
	//Array should already be hashed, sorted, unique, subtracted from Long.MAX_VALUE, then reversed.
	public Sketch(long[] array_){
		this(array_, -1, -1, -1, -1, -1, null, null, null);
	}
	
	public Sketch(SketchHeap heap, boolean clearFname){
		this(SketchTool.toSketchArray(heap), (int)heap.taxID, heap.imgID, heap.genomeSizeBases, heap.genomeSizeKmers, heap.genomeSequences, 
				heap.taxName(), heap.name0(), heap.fname());
		heap.clearSet();
		heap.clear(clearFname);
//		System.err.println("size="+size+", genome="+this.genomeSize+", m"); : (int)(2+maxGenomeFraction*heap.genomeSize)+", "+this.array.length;
//		assert(false) : (int)(2+maxGenomeFraction*heap.genomeSize)+", "+this.array.length;
	}

	public Sketch(long[] array_, int taxID_, long imgID_, long gSizeBases_, long gSizeKmers_, long gSequences_, String taxName_, String name0_, String fname_){
//		assert(taxID_>0) : name0_; //TODO: remove.  //123
		array=array_;
		taxID=taxID_;
		imgID=imgID_;
		genomeSizeBases=gSizeBases_;
		genomeSizeKmers=gSizeKmers_;
		genomeSequences=gSequences_;
		
		taxName=taxName_;
		name0=name0_;
		fname=fname_;
		
		if(ImgRecord.imgMap!=null && imgID>=0 && taxID<0){
			ImgRecord record=ImgRecord.imgMap.get(imgID);
			if(record!=null){
				if(record.name!=null && taxName==null){taxName=record.name;}
				taxID=record.taxID;
			}
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	public int countMatches(Sketch other, CompareBuffer buffer, AbstractBitSet present, boolean fillPresent){
		return countMatches(array, other.array, buffer, present, fillPresent);
	}
	
	public void add(Sketch other, int maxlen){
		final long[] a=array;
		final long[] b=other.array;
		if(maxlen<1){
			assert(false);
			maxlen=1000000;
		}
		LongList list=new LongList(Tools.min(maxlen, a.length+b.length));
		
		for(int i=0, j=0; i<a.length && j<b.length; ){
			final long ka=a[i], kb=b[j];
			if(ka==kb){//match
				list.add(ka);
				i++;
				j++;
			}else if(ka<kb){
				list.add(ka);
				i++;
			}else{
				list.add(kb);
				j++;
			}
			if(list.size()>=maxlen){break;}
		}
		
		if(array.length==list.size()){
			for(int i=0; i<list.size; i++){
				array[i]=list.array[i];
			}
		}else{
			array=list.toArray();
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Comparison          ----------------*/
	/*--------------------------------------------------------------*/
	
	public static final int countMatches(long[] a, long[] b, CompareBuffer buffer, AbstractBitSet present, boolean fillPresent){
//		assert(false) : (buffer==null)+", "+fillPresent+", "+present.cardinality();
		assert(a.length>0 && b.length>0);
		int matches=0, multiMatches=0, noHits=0;
		int contamHits=0, multiContamHits=0;
		int i=0, j=0;
		assert(present==null || present.capacity()==a.length);
//		assert(false) : buffer.rbs.capacity()+", "+buffer.rbs+", "+present;
		if(present!=null){
			if(fillPresent){
				for(; i<a.length && j<b.length; ){
					final long ka=a[i], kb=b[j];
					if(ka==kb){
						present.increment(i);
						matches++;
						i++;
						j++;
					}else if(ka<kb){
						i++;
					}else{
						j++;
					}
				}
			}else{
				for(; i<a.length && j<b.length; ){
					final long ka=a[i], kb=b[j];
					if(ka==kb){
						final int count=present.getCount(i);
						if(count>1){
							multiMatches++;
						}
						
						matches++;
						i++;
						j++;
					}else if(ka<kb){
						final int count=present.getCount(i);
						if(count>0){
							contamHits++;
							if(count>1){
								multiContamHits++;
							}
						}else{
							noHits++;
						}
						
						i++;
					}else{
						j++;
					}
				}
			}
		}else{
			for(; i<a.length && j<b.length; ){
				final long ka=a[i], kb=b[j];
				if(ka==kb){
					matches++;
					i++;
					j++;
				}else if(ka<kb){
					i++;
				}else{
					j++;
				}
			}
		}
		if(buffer!=null){
			buffer.set(matches, multiMatches, noHits, contamHits, multiContamHits, i, j, a.length, b.length);
		}
		return matches;
	}
	
//	public float identity(Sketch b, float[] ret){
//		if(ret!=null){Arrays.fill(ret, 0);}
//		return identityWeighted(array, b.array, ret);
//	}
//	
//	public static float identity(long[] a, long[] b){
//		int matches=countMatches(a, b);
//		return matches/(float)(Tools.max(1, Tools.min(a.length, b.length)));
//	}
	
	@Override
	public int hashCode(){
		long gSize=genomeSizeKmers>0 ? genomeSizeKmers : genomeSizeBases;
		int code=(int) ((gSize^taxID^imgID^(name0==null ? 0 : name0.hashCode()))&Integer.MAX_VALUE);
//		System.err.println(code+", "+gSize+", "+taxID+", "+imgID+", "+name0);
		return code;
	}
	
	@Override
	public int compareTo(Sketch b){
		if(this==b){return 0;}
		if(taxID>-1 && b.taxID>-1){return taxID-b.taxID;}
		int x=taxName.compareTo(b.taxName);
		if(x!=0){return x;}
		if(name0!=null && b.name0!=null){return name0.compareTo(b.name0);}
		return name0!=null ? 1 : b.name0!=null ? -1 : 0;
	}
	
	@Override
	public boolean equals(Object b){
		if(this==b){return true;}
		if(b==null || this.getClass()!=b.getClass()){return false;}
		return equals((Sketch)b);
	}
	
	public boolean equals(Sketch b){
		return compareTo(b)==0;
	}
	
	public ByteBuilder toHeader(){
		ByteBuilder sb=new ByteBuilder();
		return toHeader(sb);
	}
	
	public ByteBuilder toHeader(ByteBuilder sb){
		sb.append("#SZ:").append(array.length);
		sb.append("\tCD:");
		sb.append(codingArray[CODING]);
		if(delta){sb.append('D');}
		if(amino){sb.append('M');}
		if(amino8){sb.append('8');}

		if(k!=defaultK || k2!=0){
			sb.append("\tK:").append(k);
			if(k2!=0){sb.append(",").append(k2);}
		}
		if(HASH_VERSION>1){sb.append("H:").append(HASH_VERSION);}
		
		if(genomeSizeBases>0){sb.append("\tGS:").append(genomeSizeBases);}
		if(genomeSizeKmers>0){sb.append("\tGK:").append(genomeSizeKmers);}
		final long ge=genomeSizeEstimate();
		if(ge>0){sb.append("\tGE:").append(ge);}
		if(genomeSequences>0){sb.append("\tGQ:"+genomeSequences);}
		if(taxID>=0){sb.append("\tID:").append(taxID);}
		if(imgID>=0){sb.append("\tIMG:").append(imgID);}
		if(spid>0){sb.append("\tSPID:").append(spid);}
		if(fname!=null){sb.append("\tFN:").append(fname);}
		if(taxName!=null){sb.append("\tNM:").append(taxName);}
		if(name0!=null){sb.append("\tNM0:").append(name0);}
		return sb;
	}
	
	public ByteBuilder toBytes(){
		return toBytes(new ByteBuilder());
	}
	
	public ByteBuilder toBytes(ByteBuilder sb){
		long prev=0;
		toHeader(sb);
		sb.append("\n");
		byte[] temp=null;
		if(CODING==A48){temp=new byte[12];}
		for(int i=0; i<array.length; i++){
			long key=array[i];
			long x=key-prev;
			if(CODING==A48){
				appendA48(x, sb, temp);
				sb.append('\n');
			}else if(CODING==HEX){
				sb.append(Long.toHexString(x)).append('\n');
			}else if(CODING==RAW){
				sb.append(x).append('\n');
			}else{
				assert(false);
			}
			if(delta){prev=key;}
		}
		return sb;
	}
	
	public static final void appendA48(long value, ByteBuilder sb, byte[] temp){
		int i=0;
//		long value=value0;
		while(value!=0){
			byte b=(byte)(value&0x3F);
//			assert(i<temp.length) : i+", "+temp.length+", "+value0;
			temp[i]=b;
			value=value>>6;
			i++;
		}
		if(i==0){
			sb.append((byte)'0');
		}else{
			for(i--;i>=0;i--){
				sb.append((char)(temp[i]+48));
			}
		}
	}
	
	public String toString(){
		return toBytes().toString();
	}
	
	public static long parseA48(String line){
		if(line.length()==0){return 0;}
		long x=0;
		for(int i=0; i<line.length(); i++){
			x<<=6;
			x|=(((long)line.charAt(i))-48);
		}
		return x;
	}
	
	public static long parseHex(String line){
		if(line.length()==0){return 0;}
		long x=0;
		for(int i=0; i<line.length(); i++){
			x<<=4;
			x|=hexTable[line.charAt(i)];
		}
		if(line.charAt(0)=='-'){x*=-1;}
		return x;
	}
	
	public static long parseA48(byte[] line){
		if(line.length==0){return 0;}
		long x=0;
		for(byte b : line){
			x<<=6;
			x|=(((long)b)-48);
		}
		return x;
	}
	
	public static long parseHex(byte[] line){
		if(line.length==0){return 0;}
		long x=0;
		for(byte b : line){
			x<<=4;
			x|=hexTable[b];
		}
		if(line[0]=='-'){x*=-1;}
		return x;
	}
	
	public long genomeSizeEstimate() {
		return array.length==0 ? 0 : Tools.min(genomeSizeKmers, genomeSizeEstimate(array[array.length-1], array.length));
	}
	
	public String name(){return taxName!=null ? taxName : name0!=null ? name0 : fname;}
	public String taxName(){return taxName;}
	public String name0(){return name0;}
	public String fname(){return fname;}
	public int length(){return array.length;}
	public void setTaxName(String s){taxName=s;}
	public void setName0(String s){name0=s;}
	public void setFname(String s){fname=s;}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	public void makeBitSets(boolean printContam, boolean index){
		assert(compareBitSet==null && indexBitSet==null);
		if(!printContam){return;}
		compareBitSet=AbstractBitSet.make(length(), bitSetBits);
		if(index){indexBitSet=AbstractBitSet.make(length(), bitSetBits);}
	}
	
	public void addToBitSet(AbstractBitSet rbs){
		compareBitSet.add(rbs);
	}
	
	public AbstractBitSet compareBitSet(){return compareBitSet;}
	
	
	public AbstractBitSet indexBitSet(){return indexBitSet;}
	
	public void mergeBitSets(){
		assert(!mergedBitSets);
		if(compareBitSet!=null && indexBitSet!=null){
			compareBitSet.setToMax(indexBitSet);
		}
		indexBitSet=null;
		mergedBitSets=true;
	}
	
	public boolean merged(){return mergedBitSets;}
	
	public long[] array;
	public int taxID;
	public final long genomeSequences;
	public final long genomeSizeBases;
	public final long genomeSizeKmers;
	private AbstractBitSet compareBitSet; //Used for comparison
	private AbstractBitSet indexBitSet;
	private String taxName;
	private String name0;
	private String fname;

	//Extended information
	public long imgID=-1;
	public long spid=-1;
//	public String seqUnitName=null;
	
	private boolean mergedBitSets=false; //Temporary for debugging
	
	
}
