package sketch;

import java.util.Locale;

import shared.Tools;

public final class Comparison extends SketchObject implements Comparable<Comparison> {
	
	/*--------------------------------------------------------------*/
	/*----------------         Constructors         ----------------*/
	/*--------------------------------------------------------------*/
	
	public Comparison(){}
	
	public Comparison(CompareBuffer buffer){
		this(buffer, null, null);
	}
	
	public Comparison(Sketch a_, Sketch b_){
		this(null, a_, b_);
	}
	
	public Comparison(CompareBuffer buffer, Sketch a_, Sketch b_){
		
		a=a_;
		b=b_;
		
		if(buffer!=null){setFrom(buffer);}
		
		if(b!=null){
			taxName=b.taxName();
			taxID=b.taxID;
		}

//		System.err.println(this);
//		System.err.println(b.present);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Mutators           ----------------*/
	/*--------------------------------------------------------------*/
	
	public void setFrom(CompareBuffer buffer){
		hits=buffer.matches;
		multiHits=buffer.multiMatches;
		noHits=buffer.noHits;
		
		contamHits=buffer.contamHits;
		multiContamHits=buffer.multiContamHits;
		
		refDivisor=buffer.refDivisor;
		queryDivisor=buffer.queryDivisor;
		
		refSize=buffer.refSize;
		querySize=buffer.querySize;
	}
	
	public void recompare(CompareBuffer buffer){
		assert(a.merged());
		int x=a.countMatches(b, buffer, a.compareBitSet(), false);
		assert(x==hits);
		setFrom(buffer);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public int compareTo(Comparison b) {
		
		{
			float pa=score(), pb=b.score();
			if(pa>pb){
				return 1;
			}else if (pa<pb){
				return -1;
			}
		}
		
		int x=hits-b.hits;
		if(x!=0){return x;}
		x=b.minDivisor()-minDivisor();
		if(x!=0){return x;}
		x=b.maxDivisor()-maxDivisor();
		if(x!=0){return x;}
		x=b.refDivisor-refDivisor;
		if(x!=0){return x;}
		x=taxID()-b.taxID();
		if(x!=0){return x;}
		if(name0()!=null && b.name0()!=null){
			return name0().compareTo(b.name0());
		}
		if(taxName()!=null && b.taxName()!=null){
			return taxName().compareTo(b.taxName());
		}
		return 0;
	}
	
	public boolean equals(Object b){
		if(b==null || b.getClass()!=this.getClass()){return false;}
		return compareTo((Comparison)b)==0;
	}
	
	//WKID
	public float idMinDivisor(){
		return hits/(float)minDivisor();
	}
	
	//KID
	public float idMaxDivisor(){
		return hits/(float)maxDivisor();
	}
	
	public float idQueryDivisor(){
		return hits/(float)(Tools.max(1, refDivisor));
	}
	
	public float idRefDivisor(){
		return hits/(float)(Tools.max(1, refDivisor));
	}
	
	public float completeness(){
		float complt=Tools.min(1, (queryDivisor-contamHits)/(float)Tools.max(1, refDivisor));
		return complt;
//		float c2=hits/(float)refDivisor;
//		assert(queryDivisor-contamHits>=hits);
//		assert(c1>=c2);
//		System.err.println(hits+", "+contamHits+", "+refDivisor+", "+queryDivisor+", "+c1+", "+c2);
//		return Tools.max(c1, c2);
//		float kid=idMaxDivisor(), wkid=idMinDivisor();
//		return kid==0 ? 0 : kid/wkid;
	}
	
	public float contamFraction(){
		return Tools.min(1, contamHits/(float)Tools.max(1, queryDivisor));
	}
	
	public float uContamFraction() {
		int uContamHits=contamHits-multiContamHits;
		return Tools.min(1, uContamHits/(float)Tools.max(1, queryDivisor));
	}
	
	public float ani(){
		if(hits<1){return 0;}
		
		double wkid=idMinDivisor();
		return wkidToAni(wkid);

//		final float rID=hits/(float)(refDivisor);
//		final float qID=hits/(float)(queryDivisor-contamHits);
//		final float wkid2=Tools.max(qID, rID);
//		final float ani=wkidToAni(wkid2);
//		
////		System.err.println("rid: "+wkidToAni(rID)+", qid: "+wkidToAni(qID)+", qid2: "+wkidToAni(hits/(float)(queryDivisor)));
//		
//		return ani;
	}
	
	int minDivisor(){return Tools.max(1, Tools.min(refDivisor, queryDivisor));}
	int maxDivisor(){return Tools.max(1, refDivisor, queryDivisor);}
	
	public float score(){
		long est=useSizeEstimate ? genomeSizeEstimate() : genomeSizeKmers();
		float wkid=idMinDivisor();
		float kid=idMaxDivisor();
		float complt=completeness();
		float contam=contamFraction();
		return (float)Math.sqrt(80*(40000+hits+uHits())*(wkid*kid*Math.pow(est*complt, 0.2)*(1-contam*0.1)));
	}
	
	public String score2(){
		float x=score();
		if(x>999){
			return(""+(long)Math.round(x));
		}else if(x>99){
			return String.format(Locale.ROOT, "%.1f", x);
		}else if(x>9){
			return String.format(Locale.ROOT, "%.2f", x);
		}
		return String.format(Locale.ROOT, "%.3f", x);
	}
	
	public String toString(){
		return hits+", "+refDivisor+", "+queryDivisor+", "+refSize+", "+querySize+", "+contamHits;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Getters           ----------------*/
	/*--------------------------------------------------------------*/

	public String name(){return taxName!=null ? taxName : name0()!=null ? name0() : fname();}
	public String taxName(){return taxName;}
	String name0(){return b.name0();}
	String fname(){return b.fname();}

//	public int taxID(){return b.taxID<minFakeID ? b.taxID : 0;}
	public int taxID(){return (taxID<minFakeID && taxID>=0) ? taxID : 0;}
	long imgID(){return (b.imgID>0 ? b.imgID : 0);}
	
	long genomeSizeBases(){return b.genomeSizeBases;}
	long genomeSizeKmers(){return b.genomeSizeKmers;}
	long genomeSequences(){return b.genomeSequences;}
	long genomeSizeEstimate(){return b.genomeSizeEstimate();}

	public int uHits() {return hits-multiHits;}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	public Sketch a, b;

	String taxName;
	int taxID;
	
	int hits;
	int multiHits;
	int noHits;
	
	int contamHits;
	int multiContamHits;
	
	int refDivisor;
	int queryDivisor;
	
	int refSize;
	int querySize;
	
}
