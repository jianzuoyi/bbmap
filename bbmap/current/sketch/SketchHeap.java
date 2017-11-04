package sketch;

import shared.Tools;
import structures.LongHashMap;
import structures.LongHeap;
import structures.LongHeapMap;
import structures.LongHeapSet;
import structures.LongHeapSetInterface;

public class SketchHeap {
	
	SketchHeap(int limit, int minKeyOccuranceCount_){
		minKeyOccuranceCount=minKeyOccuranceCount_;
		setMode=minKeyOccuranceCount<2;
		if(setMode){
			setOrMap=set=new LongHeapSet(limit);
			map=null;
			heap=set.heap;
		}else{
			limit=(int)Tools.min(100000000, limit*(long)SketchObject.sketchHeapFactor);
			setOrMap=map=new LongHeapMap(limit);
			set=null;
			heap=map.heap;
		}
	}
	
	public void clear(boolean clearFname){
		taxID=-1;
		imgID=-1;
		genomeSizeBases=0;
		genomeSizeKmers=0;
		genomeSequences=0;
		taxName=null;
		name0=null;
		if(clearFname){fname=null;}
		setOrMap.clear();
	}
	
	public void add(SketchHeap b){
		if(taxID<0){taxID=b.taxID;}
		if(imgID<0){imgID=b.imgID;}
		if(taxName==null){taxName=b.taxName;}
		if(name0==null){name0=b.name0;}
		if(fname==null){fname=b.fname;}
		genomeSizeBases+=b.genomeSizeBases;
		genomeSizeKmers+=b.genomeSizeKmers;
		genomeSequences+=b.genomeSequences;
		if(setMode){
			set.add(b.set);
		}else{
			map.add(b.map);
		}
	}
	
	public StringBuilder toHeader(){
		StringBuilder sb=new StringBuilder();
		sb.append("#SZ:"+setOrMap.size());
		
		sb.append("\tCD:");
		sb.append(SketchObject.codingArray[SketchObject.CODING]);
		if(SketchObject.delta){sb.append('D');}
		if(SketchObject.amino){sb.append('M');}
		if(SketchObject.amino8){sb.append('8');}
		
		if(SketchObject.k!=SketchObject.defaultK){sb.append("\tK:").append(SketchObject.k);}
		if(SketchObject.HASH_VERSION>1){sb.append("H:").append(SketchObject.HASH_VERSION);}

		if(genomeSizeBases>0){sb.append("\tGS:"+genomeSizeBases);}
		if(genomeSizeKmers>0){sb.append("\tGK:"+genomeSizeKmers);}
		final long ge=genomeSizeEstimate();
		if(ge>0){sb.append("\tGE:").append(ge);}
		if(genomeSequences>0){sb.append("\tGQ:"+genomeSequences);}
		if(taxID>=0){sb.append("\tID:"+taxID);}
		if(imgID>=0){sb.append("\tIMG:"+imgID);}
		if(taxName!=null){sb.append("\tNM:"+taxName);}
		if(name0!=null){sb.append("\tNM0:"+name0);}
		if(fname!=null){sb.append("\tFN:"+fname);}
		return sb;
	}
	
	public boolean checkAndAdd(long value){
		assert(value>=SketchObject.minHashValue);
		
//		if(!heap.hasRoom() && value<=heap.peek()){return false;}
//		if(Blacklist.contains(value)){return false;}
//		if(!Whitelist.contains(value)){return false;}
		
		if(Blacklist.exists() || Whitelist.exists()){
			if(!heap.hasRoom() && value<=heap.peek()){return false;}
			if(Blacklist.contains(value)){return false;}
			if(!Whitelist.containsRaw(value)){return false;}
		}
		
		return add(value);
	}
	
	@Override
	public int hashCode(){
		long gSize=genomeSizeKmers>0 ? genomeSizeKmers : genomeSizeBases;
		int code=(int) ((gSize^taxID^imgID^(name0==null ? 0 : name0.hashCode()))&Integer.MAX_VALUE);
		return code;
	}
	
	public long genomeSizeEstimate() {
		int size=size();
		if(size==0){return 0;}
		long min=peek();
		return Tools.min(genomeSizeKmers, SketchObject.genomeSizeEstimate(Long.MAX_VALUE-min, size));
	}
	
	public long sketchSizeEstimate(){
		return SketchObject.toSketchSize(genomeSizeBases, genomeSizeKmers, genomeSizeEstimate(), SketchObject.targetSketchSize);
	}
	
	public String toString(){return toHeader().toString();}

	public String name(){return taxName==null ? name0 : taxName;}
	public String taxName(){return taxName;}
	public String name0(){return name0;}
	public String fname() {return fname;}
	public void setTaxName(String s){taxName=s;}
	public void setName0(String s){name0=s;}
	public void setFname(String s) {fname=s;}
	
	private String taxName;
	private String name0;
	private String fname;
	public long taxID=-1;
	public long imgID=-1;
	public long genomeSizeBases=0;
	public long genomeSizeKmers=0;
	public long genomeSequences=0;

	public int capacity(){return heap.capacity();}
	public boolean hasRoom(){return heap.hasRoom();}
	public long peek(){return heap.peek();}
	public int size(){return heap.size();}
	public LongHashMap map(){return map.map;}

	public void clear(){setOrMap.clear();}
	public void clearSet(){
		if(set==null){map.map.clear();}
		else{set.set.clear();}
	}
	public boolean add(long key){return setOrMap.add(key);}

	private final LongHeapSet set;
	private final LongHeapMap map;
	private final LongHeapSetInterface setOrMap;
	public final LongHeap heap;
	public final int minKeyOccuranceCount;
	public final boolean setMode;
}
