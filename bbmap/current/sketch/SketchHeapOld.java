package sketch;

import shared.Tools;
import structures.LongHeapSet;

public class SketchHeapOld extends LongHeapSet {
	
	SketchHeapOld(int limit_){
		super(limit_);
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
		super.clear();
	}
	
	public void add(SketchHeapOld b){
		if(taxID<0){taxID=b.taxID;}
		if(imgID<0){imgID=b.imgID;}
		if(taxName==null){taxName=b.taxName;}
		if(name0==null){name0=b.name0;}
		if(fname==null){fname=b.fname;}
		genomeSizeBases+=b.genomeSizeBases;
		genomeSizeKmers+=b.genomeSizeKmers;
		genomeSequences+=b.genomeSequences;
		super.add(b);
	}
	
	public StringBuilder toHeader(){
		StringBuilder sb=new StringBuilder();
		sb.append("#SZ:"+size());
		
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
		
		return super.add(value);
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
		long min=heap.peek();
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
	
}
