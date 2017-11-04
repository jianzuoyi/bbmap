package sketch;

import shared.Tools;
import structures.AbstractBitSet;

public class CompareBuffer {
	
	public CompareBuffer(boolean makeBS){
		if(makeBS){
			cbs=AbstractBitSet.make(0, SketchObject.bitSetBits);
		}else{
			cbs=null;
		}
	}
	
	void set(final int matches_, final int multiMatches_, final int noHits_, final int contamHits_, final int multiContamHits_,
			final int queryDivisor_, final int refDivisor_, final int querySize_, final int refSize_){
		matches=matches_;
		multiMatches=multiMatches_;
		noHits=noHits_;
		
		contamHits=contamHits_;
		multiContamHits=multiContamHits_;
		
		queryDivisor=queryDivisor_;
		refDivisor=refDivisor_;
		
		querySize=querySize_;
		refSize=refSize_;
	}

	void clear(){
		matches=multiMatches=noHits=0;
		contamHits=multiContamHits=0;
		refDivisor=queryDivisor=0;
		refSize=querySize=0;
	}

	int minDivisor(){return Tools.max(1, Tools.min(queryDivisor, refDivisor));}
	int maxDivisor(){return Tools.max(1, queryDivisor, refDivisor);}
	int minSize(){return Tools.max(1, Tools.min(querySize, refSize));}
	int maxSize(){return Tools.max(1, querySize, refSize);}

	int uniqueMatches(){return matches-multiMatches;}
	int uniqueContamHits(){return contamHits-multiContamHits;}
	
	int matches;
	int multiMatches;
	int noHits;
	
	int contamHits;
	int multiContamHits;
	
	int queryDivisor;
	int refDivisor;
	
	int querySize;
	int refSize;

	public final AbstractBitSet cbs; //Only for comparisons, not index
	
}
