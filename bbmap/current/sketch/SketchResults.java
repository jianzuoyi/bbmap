package sketch;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Map.Entry;
import java.util.concurrent.ConcurrentHashMap;

import fileIO.TextStreamWriter;
import shared.Shared;
import structures.IntHashMap;

public class SketchResults extends SketchObject {
	
	SketchResults(Sketch s){
		sketch=s;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	public void addMap(ConcurrentHashMap<Integer, Comparison> map, DisplayParams params, CompareBuffer buffer) {

		if(map.isEmpty()){return;}
		list=addToList(map, params, list);
		
		if((params.printContam || params.printContamHits || true)/* && (makeIndex || !AUTOSIZE)*/){
			recompare(buffer);
		}
	}
	
	public void recompare(CompareBuffer buffer){
//		assert(makeIndex || !AUTOSIZE);
		
		assert(!sketch.merged());
		sketch.mergeBitSets();
		
//		System.err.println(sketch.compareBitSet());
//		assert(false) : sketch.compareBitSet().getClass();
		
		for(Comparison c : list){
			c.recompare(buffer);
		}
		Collections.sort(list);
		Collections.reverse(list);
	}
	
	private static ArrayList<Comparison> addToList(ConcurrentHashMap<Integer, Comparison> map, DisplayParams params, ArrayList<Comparison> old){

//		System.err.println(map.size());
//		System.err.println(map.keySet());
		
		ArrayList<Comparison> al=(old==null ? new ArrayList<Comparison>(map.size()) : old);
		for(Entry<Integer, Comparison> e : map.entrySet()){
			al.add(e.getValue());
		}
		Shared.sort(al);
		Collections.reverse(al);
		while(al.size()>params.maxRecords*2+10){
			al.remove(al.size()-1);
		}
		return al;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Tax Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	public int primaryTax(int level){
		//I have no idea how to implement this...
		IntHashMap map=new IntHashMap();
		assert(false);
		return -1;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Print Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	private static String recordBreak="\n"; //"\n\n"
	
	void writeResults(DisplayParams params, TextStreamWriter tsw){
		StringBuilder sb=toText(params);
		tsw.print(sb);
	}

	StringBuilder toText(DisplayParams params){
		final StringBuilder sb=params.queryHeader(sketch);
		if(params.format==3){
			if(list==null || list.isEmpty()){return sb;}
			int idx=0;
			int prevTaxID=0;
			for(Comparison c : list){
				params.formatComparison(c, sb, prevTaxID, sketch);
				prevTaxID=c.taxID();
				idx++;
				if(idx>=params.maxRecords){break;}
			}
		}else{
			sb.append(recordBreak);

			if(list==null || list.isEmpty()){
				sb.append("No hits.\n");
			}else{
				if(params.format==2){sb.append(params.header()).append('\n');}
				int idx=0;
				int prevTaxID=0;
				for(Comparison c : list){
					params.formatComparison(c, sb, prevTaxID, sketch);
					prevTaxID=c.taxID();
					idx++;
					if(idx>=params.maxRecords){break;}
				}
			}
		}
		return sb;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	public final Sketch sketch;
	public ArrayList<Comparison> list;
	public int totalRecords=0;
	
}
