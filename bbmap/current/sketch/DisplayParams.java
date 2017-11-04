package sketch;

import java.util.ArrayList;
import java.util.Locale;

import shared.Colors;
import shared.Tools;
import tax.PrintTaxonomy;
import tax.TaxNode;
import tax.TaxTree;

public class DisplayParams implements Cloneable {
	
	@Override
	public DisplayParams clone(){
		try {
			return (DisplayParams)super.clone();
		} catch (CloneNotSupportedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			throw new RuntimeException();
		}
	}
	
	public DisplayParams parseDoubleHeader(String s){
		if(!s.startsWith("##")){return this;}
		StringBuilder sb=new StringBuilder();
		for(int i=2; i<s.length(); i++){
			char c=s.charAt(i);
			if(c=='\n'){break;}
			sb.append(c);
		}
		return parseDoubleHeaderLine(sb.toString());
	}
	
	public DisplayParams parseDoubleHeaderLine(String line) {
		if(line.startsWith("##")){line=line.substring(2);}
		else{assert(!line.startsWith("#")) : line;}
		if(line.length()<1){return this;}
		
		DisplayParams params=(DisplayParams) this.clone();
		
		String[] args=line.split(" ");
		for(String arg : args){
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b==null || b.equalsIgnoreCase("null")){b=null;}
			while(a.startsWith("-")){a=a.substring(1);} //Strip leading hyphens
			
			boolean x=params.parse(arg, a, b);
//			assert(x) : "Unknown parameter "+arg+"\n"+line;
			if(!x){System.err.println("Warning: Unknown parameter "+arg);}
		}
		
		return params;
	}
	
	public boolean parse(String arg, String a, String b){
	
		if(a.equals("minhits")  || a.equals("hits")){
			minHits=Integer.parseInt(b);
		}else if(a.equalsIgnoreCase("minwkid") || a.equalsIgnoreCase("wkid")){
			minWKID=Float.parseFloat(b);
			if(minWKID>1){minWKID/=100;}
			assert(minWKID<=1) : "minWKID should between 0 and 1";
		}else if(a.equalsIgnoreCase("minid") || a.equalsIgnoreCase("id") || a.equalsIgnoreCase("minani") || a.equalsIgnoreCase("ani")){
			minANI=Float.parseFloat(b);
			if(minANI>1){minANI/=100;}
			assert(minANI<=1) : "minANI should between 0 and 1";
			if(minANI>0){
				minWKID=(float)Tools.max(minWKID, Comparison.aniToWkid(minANI, 32));//Lowest possible minWKID for this ANI
			}
		}else if(a.equals("records") || a.equals("maxrecords") || a.equals("results")){
			maxRecords=Integer.parseInt(b);
			assert(maxRecords>=1) : "Max records must be at least 1.";
		}else if(a.equals("format")){
			format=Integer.parseInt(b);
		}else if(a.equals("level") || a.equals("taxlevel") || a.equals("minlevel")){
			if(b==null){taxLevel=-1;}
			else if(Character.isDigit(b.charAt(0))){
				taxLevel=Integer.parseInt(b);
			}else{
				taxLevel=TaxTree.stringToLevel(b.toLowerCase());
			}
		}
		
		else if(a.equalsIgnoreCase("printtax") || a.equalsIgnoreCase("printtaxa")){
			printTax=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printoriginalname") || a.equalsIgnoreCase("printseqname") || a.equalsIgnoreCase("printname0") || a.equals("pn0")){
			printOriginalName=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printfilename") || a.equalsIgnoreCase("printfname")){
			printFileName=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printimg")){
			printImg=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printcompleteness") || a.equalsIgnoreCase("completeness") || a.equalsIgnoreCase("printcomplt")){
			printCompleteness=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printani") || a.equalsIgnoreCase("ani")){
			printAni=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printscore") || a.equalsIgnoreCase("score")){
			printScore=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("contamhits") || a.equalsIgnoreCase("contam") || a.equalsIgnoreCase("printcontam")){
			printContam=Tools.parseBoolean(b);
		}
		
		else if(a.equalsIgnoreCase("printMatches")){
			printMatches=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printLength")){
			printLength=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printTaxID")){
			printTaxID=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printGSize")){
			printGSize=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printGKmers")){
			printGKmers=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printTaxName")){
			printTaxName=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printGSeqs")){
			printGSeqs=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printGBases")){
			printGBases=Tools.parseBoolean(b);
		}
		
		else if(a.equalsIgnoreCase("printUMatches") || a.equalsIgnoreCase("printUHits") || a.equalsIgnoreCase("printUnique")){
			printUMatches=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printUContam")){
			printUContam=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printNoHit")){
			printNoHit=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printColors") || a.equalsIgnoreCase("colors") || a.equalsIgnoreCase("color")){
//			System.err.println("Parsing '"+arg+"'"); //123
			if(b==null || b.length()<1){
				printColors=true;
			}else if(b.equalsIgnoreCase("t") || b.equalsIgnoreCase("true")){
				printColors=true;
			}else if(b.equalsIgnoreCase("f") || b.equalsIgnoreCase("false")){
				printColors=false;
			}else{
				printColors=true;
				if(Character.isDigit(b.charAt(0)) || b.charAt(0)=='-'){
					colorLevel=Tools.max(0, TaxTree.levelToExtended(Integer.parseInt(b)));
				}else{
					colorLevel=TaxTree.stringToLevelExtended(b);
				}
			}
			setColors=true;
//			System.err.println("Parsed "+arg); //123
		}else if(a.equalsIgnoreCase("colorLevel")){
//			System.err.println("Parsing '"+arg+"'"); //123
			if(Character.isDigit(b.charAt(0)) || b.charAt(0)=='-'){
				colorLevel=Tools.max(0, TaxTree.levelToExtended(Integer.parseInt(b)));
			}else{
				colorLevel=TaxTree.stringToLevelExtended(b);
			}
//			System.err.println("Parsed "+arg); //123
		}
		
		else if(a.equalsIgnoreCase("printRefDivisor") || a.equalsIgnoreCase("printRDiv")){
			printRefDivisor=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printQueryDivisor") || a.equalsIgnoreCase("printQDiv")){
			printQueryDivisor=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printRefSize") || a.equalsIgnoreCase("printRSize")){
			printRefSize=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printQuerySize") || a.equalsIgnoreCase("printQSize")){
			printQuerySize=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printContamHits") || a.equalsIgnoreCase("printCHits")){
			printContamHits=Tools.parseBoolean(b);
		}
		
		else if(a.equalsIgnoreCase("printAll")){
			if(Tools.parseBoolean(b)){
				setPrintAll();
			}
		}
		
		else if(a.equals("samplerate")){
			samplerate=Float.parseFloat(b);
		}else if(a.equals("reads")){
			reads=Tools.parseKMG(b);
		}else if(a.equals("mode")){
			mode=SketchObject.parseMode(arg, a, b);
		}
		
		//For format 3
		else if(a.equalsIgnoreCase("useTaxidName") || a.equalsIgnoreCase("useTaxidAsName")){
			useTaxidName=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("useImgName") || a.equalsIgnoreCase("useImgAsName")){
			useImgName=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("useTaxName") || a.equalsIgnoreCase("useTaxAsName")){
			useTaxName=Tools.parseBoolean(b);
		}
		
		else if(a.equalsIgnoreCase("minkmercount") || a.equalsIgnoreCase("minkeycount") || a.equalsIgnoreCase("minKeyOccuranceCount")){
			minKeyOccuranceCount=Integer.parseInt(b);
		}
		
		//TODO:  Eventually remove support for "amino" and "k" and just support "hamino" and "hk"
		
		//Parameters for compatibility verification
		else if(a.equalsIgnoreCase("k") || a.equalsIgnoreCase("hk")){
//			System.err.println("A: k="+k+", k2="+k2+", arg="+arg);
			if(b.indexOf(',')>=0){
				String[] split=b.split(",");
				assert(split.length==2) : "\nBad argument "+arg+"\n"+b+"\n";
				int x=Integer.parseInt(split[0]);
				int y=Integer.parseInt(split[1]);
				k=Tools.max(x, y);
				k2=Tools.min(x, y);
//				System.err.println("B: k="+k+", k2="+k2+", split="+Arrays.toString(split));
			}else{
				k=Integer.parseInt(b);
//				System.err.println("C: k="+k+", k2="+k2);
			}
		}else if(a.equals("hashversion") || a.equals("hv")){
			hashVersion=Integer.parseInt(b);
		}else if(a.equals("amino") || a.equals("hamino")){
			amino=Tools.parseBoolean(b);
		}
		
		else{
			return false;
		}
		return true;
	}
	
	public String toString(){
		StringBuilder sb=new StringBuilder();
		sb.append("##");
		sb.append("hits=").append(minHits);
		sb.append(" wkid=").append(String.format("%.5f",minWKID));
		if(minANI>0){sb.append(" id=").append(String.format("%.5f",minANI));}
		sb.append(" records=").append(maxRecords);
		sb.append(" format=").append(format);
		sb.append(" level=").append(taxLevel);
		
		if(k!=SketchObject.defaultK || k2!=0 || k!=SketchObject.k || k2!=SketchObject.k2){
			assert(k>0 && k2>=0 && k2<k) : "Bad values for k: "+k+", "+k2+", "+SketchObject.k+", "+SketchObject.k2;
			assert(SketchObject.k>0 && SketchObject.k2>=0 && SketchObject.k2<SketchObject.k) : "Bad values for k: "+k+", "+k2+", "+SketchObject.k+", "+SketchObject.k2;
			sb.append(" hk=").append(SketchObject.k).append(',').append(SketchObject.k2);
		}
		if(SketchObject.amino){sb.append(" hamino=").append(SketchObject.amino);} //TODO: This conflicts with Parser flag
		if(SketchObject.HASH_VERSION>1){sb.append(" hashversion=").append(SketchObject.HASH_VERSION);}

		if(true || printTax!=default_printTax){sb.append(" printTax=").append(printTax);}
		if(true || printOriginalName!=default_printOriginalName){sb.append(" pn0=").append(printOriginalName);}
		if(true || printFileName!=default_printFileName){sb.append(" printfname=").append(printFileName);}
		if(true || printImg!=default_printImg){sb.append(" printImg=").append(printImg);}
		if(true || printAni!=default_printAni){sb.append(" printAni=").append(printAni);}
		if(true || printCompleteness!=default_printCompleteness){sb.append(" printCompleteness=").append(printCompleteness);}
		if(true || printScore!=default_printScore){sb.append(" printScore=").append(printScore);}
		if(true || printContam!=default_printContam){sb.append(" contam=").append(printContam);}
		
		if(true || printMatches!=default_printMatches){sb.append(" printMatches=").append(printMatches);}
		if(true || printLength!=default_printLength){sb.append(" printLength=").append(printLength);}
		if(true || printTaxID!=default_printTaxID){sb.append(" printTaxID=").append(printTaxID);}
		if(true || printGSize!=default_printGSize){sb.append(" printGSize=").append(printGSize);}
		if(true || printGKmers!=default_printGKmers){sb.append(" printGKmers=").append(printGKmers);}
		if(true || printTaxName!=default_printTaxName){sb.append(" printTaxName=").append(printTaxName);}
		if(true || printGSeqs!=default_printGSeqs){sb.append(" printGSeqs=").append(printGSeqs);}
		if(true || printGBases!=default_printGBases){sb.append(" printGBases=").append(printGBases);}
		
		if(true || printUMatches!=default_printUMatches){sb.append(" printUMatches=").append(printUMatches);}
		if(true || printUContam!=default_printUContam){sb.append(" printUContam=").append(printUContam);}
		if(true || printNoHit!=default_printNoHit){sb.append(" printNoHit=").append(printNoHit);}

		if(useTaxidName){sb.append(" useTaxidName=").append(useTaxidName);}
		if(useImgName){sb.append(" useImgName=").append(useImgName);}
		if(useTaxName){sb.append(" useTaxName=").append(useTaxName);}
		
		if(true){sb.append(" colors=").append(printColors ? TaxTree.extendedToLevel(colorLevel)+"" : "f");}
		
		if(minKeyOccuranceCount!=default_minKeyOccuranceCount){sb.append(" minKeyOccuranceCount=").append(minKeyOccuranceCount);}
		
//		if(printColors && colorLevel!=default_colorLevel){sb.append(" colorLevel=").append(TaxTree.extendedToLevel(colorLevel));}
		

		if(printRefDivisor){sb.append(" printRefDivisor=").append(printRefDivisor);}
		if(printQueryDivisor){sb.append(" printQueryDivisor=").append(printQueryDivisor);}
		if(printRefSize){sb.append(" printRefSize=").append(printRefSize);}
		if(printQuerySize){sb.append(" printQuerySize=").append(printQuerySize);}
		if(printContamHits){sb.append(" printContamHits=").append(printContamHits);}
		
		
		if(reads>-1){sb.append(" reads=").append(reads);}
		if(mode!=default_mode){sb.append(" mode=").append(mode);}
		if(samplerate!=default_samplerate){sb.append(" samplerate=").append(String.format("%.4f",samplerate));}
		
		
		sb.append('\n');
		return sb.toString();
	}
	
	public boolean compatible(){
		return SketchObject.k==k && SketchObject.k2==k2 && SketchObject.amino==amino && hashVersion==SketchObject.HASH_VERSION;
	}
	
	public void setPrintAll(){
		printTax=true;
		printOriginalName=true;
		printImg=true;
		printAni=true;
		printCompleteness=true;
		printScore=true;
		
		printMatches=true;
		printLength=true;
		printTaxID=true;
		printGSize=true;
		printGKmers=true;
		printTaxName=true;
		printGSeqs=true;
		printGBases=true;

		printUMatches=true;
		printUContam=true;
		printNoHit=true;
		
//		printColors=true;
		
		printContam=true;
		
		printRefDivisor=true;
		printQueryDivisor=true;
		printRefSize=true;
		printQuerySize=true;
		printContamHits=true;
	}
	

	
	/*--------------------------------------------------------------*/
	/*----------------          Formatting          ----------------*/
	/*--------------------------------------------------------------*/

	StringBuilder queryHeader(Sketch s){
		StringBuilder sb=new StringBuilder();
		if(format>2){return sb;}
		
		String color=toColor(s.taxID);
		if(color!=null){sb.append(color);}
		
		sb.append("\nQuery: ").append(s.name()==null ? "." : s.name());
		if(dbName!=null){sb.append("\tDB: ").append(dbName);}
		sb.append("\tSeqs: ").append(s.genomeSequences).append(' ');
		sb.append("\tBases: ").append(s.genomeSizeBases);
		sb.append("\tgSize: ").append(s.genomeSizeEstimate());
		sb.append("\tSketchLen: ").append(s.length());

		if(s.imgID>0){sb.append("\tIMG: ").append(s.imgID);}
		if(s.spid>0){sb.append("\tspid: ").append(s.spid);}
		if(s.taxID>0 && s.taxID<SketchObject.minFakeID){sb.append("\tTaxID: ").append(s.taxID);}

		if(printFileName && s.fname()!=null && !s.fname().equals(s.name())){sb.append("\tFile: "+s.fname());}
		if(printOriginalName && s.name0()!=null && !s.name0().equals(s.name())){sb.append("\tSeqName: "+s.name0());}
		
		if(color!=null){sb.append(Colors.RESET);}
		
		return sb;
	}
	
	int toColorTid(final int taxID){
		if(!printColors || SketchObject.taxtree==null || taxID<=0 || taxID>=SketchObject.minFakeID){return 0;}
		TaxNode tn=SketchObject.taxtree.getNode(taxID);
		while(tn!=null && tn.id!=tn.pid && tn.levelExtended<colorLevel){
			tn=SketchObject.taxtree.getNode(tn.pid);
//			System.err.println(tn);
		}
		return tn==null || tn.levelExtended>=TaxTree.LIFE_E || (tn.levelExtended>colorLevel && tn.levelExtended>TaxTree.PHYLUM_E) ? 0 : tn.id;
	}
	
	String toColor(final int taxID){
		if(!printColors || SketchObject.taxtree==null || taxID<=0 || taxID>=SketchObject.minFakeID){return null;}
		TaxNode tn=SketchObject.taxtree.getNode(taxID);
		while(tn!=null && tn.id!=tn.pid && tn.levelExtended<colorLevel){
			tn=SketchObject.taxtree.getNode(tn.pid);
//			System.err.println(tn);
		}
		if(tn==null){
			return null;
		}else{
			if(tn.levelExtended>=TaxTree.LIFE_E || (tn.levelExtended>colorLevel && tn.levelExtended>TaxTree.PHYLUM_E)){return Colors.WHITE;}
			else{
//				System.err.println("*"+tn.id+", "+tn.id%Colors.colorArray.length);
				return Colors.colorArray[tn.id%Colors.colorArray.length];
			}
		}
	}
	
	String header(){
		if(format==3){return "Query\tReference\tANI";}
		
		StringBuilder sb=new StringBuilder();
		
		//Numeric fields
		if(true){sb.append("WKID");}
		if(true){sb.append("\tKID");}
		if(printAni){sb.append("\tANI");}
		if(printCompleteness){sb.append("\tComplt");}
		if(printContam){sb.append("\tContam");}
		if(printUContam){sb.append("\tuContam");}
		if(printScore){sb.append("\tScore");}
		if(printMatches){sb.append("\tMatches");}
		if(printUMatches){sb.append("\tUnique");}
		if(printNoHit){sb.append("\tnoHit");}
		if(printLength){sb.append("\tLength");}
		if(printTaxID){sb.append("\tTaxID");}
		if(printImg){sb.append("\tImgID");}
		if(printGBases){sb.append("\tgBases");}
		if(printGKmers){sb.append("\tgKmers");}
		if(printGSize){sb.append("\tgSize");}
		if(printGSeqs){sb.append("\tgSeqs");}
		
		
		//Raw fields
		if(printRefDivisor){sb.append("\trDiv");}
		if(printQueryDivisor){sb.append("\tqDiv");}
		if(printRefSize){sb.append("\trSize");}
		if(printQuerySize){sb.append("\tqSize");}
		if(printContamHits){sb.append("\tcHits");}
		
		//Text fields
		if(printTaxName){sb.append("\ttaxName");}
		if(printOriginalName){sb.append("\tseqName");}
		if(printTax && SketchObject.taxtree!=null){sb.append("\ttaxonomy");}
		
		return sb.toString();
	}
	
	void formatComparisonColumnwise(Comparison c, StringBuilder sb, int prevTid){
		final int tid=c.taxID;
		boolean reset=false;
		
		if(printColors){
			final int ctid=toColorTid(tid);
			final int prevCtid=toColorTid(prevTid);

			final int cnum=ctid%Colors.colorArray.length;
			final int prevCnum=prevCtid%Colors.colorArray.length;

			String color=toColor(tid);
			String underline=(printColors && cnum==prevCnum && ctid!=prevCtid && (ctid>1 && prevCtid>1) ? Colors.UNDERLINE : null);

			if(color!=null){sb.append(color);}
			if(underline!=null){sb.append(underline);}
			reset=(color!=null || underline!=null);
			
//			System.err.println((color==null ? "" : color)+(underline==null ? "" : underline)+
//					tid+", "+prevTid+";     \t"+ctid+", "+prevCtid+";     \t"+cnum+", "+prevCnum+"; \t"+((underline!=null)+"")+Colors.RESET);
//			System.err.println(color==null ? "null" : color.substring(1));
		}
		
		sb.append(String.format(Locale.ROOT, "%.2f%%\t%.2f%%", 100*c.idMinDivisor(), 100*c.idMaxDivisor()));
		
		if(printAni){sb.append(String.format(Locale.ROOT, "\t%.2f%%", 100*c.ani()));}
		if(printCompleteness){sb.append(String.format(Locale.ROOT, "\t%.2f%%", 100*c.completeness()));}
		if(printContam){sb.append(String.format(Locale.ROOT, "\t%.2f%%", 100*c.contamFraction()));}
		if(printUContam){sb.append(String.format(Locale.ROOT, "\t%.2f%%", 100*c.uContamFraction()));}
		if(printScore){sb.append('\t').append(c.score2());}
		if(printMatches){sb.append('\t').append(c.hits);}
		if(printUMatches){sb.append('\t').append(c.uHits());}
		if(printNoHit){sb.append('\t').append(c.noHits);}
		if(printLength){sb.append('\t').append( c.maxDivisor());}
		if(printTaxID){sb.append('\t').append(tid);}
		if(printImg){sb.append(String.format(Locale.ROOT, "\t%d", c.imgID()));}
		if(printGBases){sb.append('\t').append(c.genomeSizeBases());}
		if(printGKmers){sb.append('\t').append(c.genomeSizeKmers());}
		if(printGSize){sb.append('\t').append(c.genomeSizeEstimate());}
		if(printGSeqs){sb.append('\t').append(c.genomeSequences());}
		
		//Raw fields
		if(printRefDivisor){sb.append('\t').append(c.refDivisor);}
		if(printQueryDivisor){sb.append('\t').append(c.queryDivisor);}
		if(printRefSize){sb.append('\t').append(c.refSize);}
		if(printQuerySize){sb.append('\t').append(c.querySize);}
		if(printContamHits){sb.append('\t').append(c.contamHits);}
		
		//Text fields
		if(printTaxName){sb.append('\t').append(c.taxName()==null ? "." : c.taxName());}
		if(printOriginalName || (c.taxName()==null && c.name0()!=null)){sb.append('\t').append(c.name0()==null ? "." : c.name0());}
		if(printTax && SketchObject.taxtree!=null){
			sb.append('\t');
			TaxNode tn=null;
			if(tid>0 && tid<SketchObject.minFakeID){
				tn=SketchObject.taxtree.getNode(tid);
			}

			if(tn!=null){
				sb.append(SketchObject.taxtree.toSemicolon(tn, SketchObject.skipNonCanonical));
			}else{
				sb.append('.');
			}
		}
		
		if(reset){sb.append(Colors.RESET);}
		
		sb.append('\n');
	}
	
	void formatComparison3Column(Comparison c, StringBuilder sb, int prevTid, Sketch query){

		final String qName=useTaxidName ? ""+query.taxID : useImgName ?  ""+query.imgID : useTaxName ? ""+query.taxName() : query.name();
		final String rName=useTaxidName ? ""+c.taxID() : useImgName ?  ""+c.imgID() : useTaxName ? ""+c.taxName() : c.name();
		final int tid=c.taxID;
		boolean reset=false;
		
		sb.append(qName).append('\t');
		
		if(printColors){
			final int ctid=toColorTid(tid);
			final int prevCtid=toColorTid(prevTid);

			final int cnum=ctid%Colors.colorArray.length;
			final int prevCnum=prevCtid%Colors.colorArray.length;

			String color=toColor(tid);
			String underline=(printColors && cnum==prevCnum && ctid!=prevCtid && (ctid>1 && prevCtid>1) ? Colors.UNDERLINE : null);

			if(color!=null){sb.append(color);}
			if(underline!=null){sb.append(underline);}
			reset=(color!=null || underline!=null);
			
//			System.err.println((color==null ? "" : color)+(underline==null ? "" : underline)+
//					tid+", "+prevTid+";     \t"+ctid+", "+prevCtid+";     \t"+cnum+", "+prevCnum+"; \t"+((underline!=null)+"")+Colors.RESET);
//			System.err.println(color==null ? "null" : color.substring(1));
		}
		
		sb.append(rName).append('\t').append(String.format(Locale.ROOT, "\t%.2f", 100*c.ani()));
		
		if(reset){sb.append(Colors.RESET);}
		
		sb.append('\n');
	}
	
	void formatComparison(Comparison c, StringBuilder sb, int prevTaxID, Sketch query){
		if(format==2){
			formatComparisonColumnwise(c, sb, prevTaxID);
			return;
		}else if(format==3){
			formatComparison3Column(c, sb, prevTaxID, query);
			return;
		}
		String complt=(printCompleteness ? String.format(Locale.ROOT, "\tcomplt %.2f%%%%", 100*c.completeness()) : "");
		String contam=(printContam ? String.format(Locale.ROOT, "\tcontam %.2f%%%%", 100*c.contamFraction()) : "");
//		String score=(printScore ? String.format(Locale.ROOT, "\tscore %.2f", c.score2()) : "");
		String score=(printScore ? "\tscore "+c.score2() : "");
		String ccs=complt+contam+score;
		
		if(format==0){
			sb.append(String.format(Locale.ROOT, "WKID %.2f%%\tKID %.2f%%"+ccs+"\tmatches %d\tcompared %d",
					100*c.idMinDivisor(), 100*c.idMaxDivisor(), c.hits, c.minDivisor())+"\ttaxID "+c.taxID()+
					(printImg ? "\timgID "+c.imgID() : "")+"\tgKmers "+c.genomeSizeKmers()+"\t"+
					(c.taxName()==null ? "." : c.taxName())+
					((printOriginalName || (c.taxName()==null && c.name0()!=null)) ? "\t"+(c.name0()==null ? "." : c.name0()) : "")+"\n");
			if(printTax && SketchObject.taxtree!=null){
				if(c.taxID()>=0 && c.taxID()<SketchObject.minFakeID){
					TaxNode tn=SketchObject.taxtree.getNode(c.taxID());
					if(tn!=null){
						PrintTaxonomy.printTaxonomy(tn, sb, SketchObject.taxtree, TaxTree.DOMAIN, SketchObject.skipNonCanonical);
					}
				}
				sb.append('\n');
			}
		}else{
			ArrayList<TaxNode> tnl=new ArrayList<TaxNode>();
			if(SketchObject.taxtree!=null && c.taxID()>=0 && c.taxID()<SketchObject.minFakeID){
				TaxNode tn=SketchObject.taxtree.getNode(c.taxID());
				while(tn!=null && tn.pid!=tn.id && tn.level<=TaxTree.DOMAIN){
					tnl.add(tn);
					tn=SketchObject.taxtree.getNode(tn.pid);
				}
			}
			
			sb.append(String.format(Locale.ROOT, "WKID %.2f%%\tKID %.2f%%"+ccs+"\tmatches %d\tcompared %d\t",
					100*c.idMinDivisor(), 100*c.idMaxDivisor(), c.hits, c.minDivisor()));
			sb.append("\ttaxID ").append(c.taxID()).append('\t');
			if(printImg){sb.append("\timgID ").append(c.imgID()).append('\t');}
			sb.append(c.taxName()).append('\t');
			if(printOriginalName || (c.taxName()==null && c.name0()!=null)){sb.append(c.name0()).append('\t');}
			
			if(printTax){
				for(int i=tnl.size()-1; i>=0; i--){
					TaxNode tn=tnl.get(i);
					sb.append(tn.name);
					if(i>0){sb.append(';');}
				}
			}
			sb.append('\n');
			
			tnl.clear();
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	//These are shared with SketchObject
	//They do not affect anything and are just for the server to validate remote settings.
	private int hashVersion=SketchObject.HASH_VERSION;
	private int k=SketchObject.k;
	private int k2=SketchObject.k2;
	private boolean amino=SketchObject.amino;
	
	boolean amino(){return amino;}
	
	//These are unique
	public int maxRecords=default_maxRecords;
	public float minANI=0;
	public float minWKID=default_minWKID;
	public int format=default_format;
	
	public int minHits=default_minHits;
	public int taxLevel=default_taxLevel;
	public int mode=default_mode;
	public float samplerate=default_samplerate;
	public long reads=default_reads;
	public int minKeyOccuranceCount=default_minKeyOccuranceCount;
	
	public String dbName=null;
	
	/*--------------------------------------------------------------*/
	/*----------------         Print Columns        ----------------*/
	/*--------------------------------------------------------------*/
	
	//For format 2
	public boolean printTax=default_printTax;
	public boolean printOriginalName=default_printOriginalName;
	public boolean printFileName=default_printFileName;
	public boolean printImg=default_printImg;
	public boolean printAni=default_printAni;
	public boolean printCompleteness=default_printCompleteness;
	public boolean printScore=default_printScore;
	
	public boolean printLength=default_printLength;
	public boolean printTaxID=default_printTaxID;
	public boolean printGSize=default_printGSize;
	public boolean printGKmers=default_printGKmers;
	public boolean printTaxName=default_printTaxName;
	public boolean printGSeqs=default_printGSeqs;
	public boolean printGBases=default_printGBases;
	
	public boolean printUMatches=default_printUMatches;
	public boolean printUContam=default_printUContam;
	public boolean printNoHit=default_printNoHit;

	public boolean printColors=default_printColors;
	public boolean setColors=false;
	public int colorLevel=default_colorLevel;
	
	/** TODO: Note this is conflated between printing %contam and calculating things based on contam hits. */
	public boolean printContam=default_printContam;
	
	/** Raw fields */
	public boolean printMatches=default_printMatches;
	
	public boolean printRefDivisor=false;
	public boolean printQueryDivisor=false;
	public boolean printRefSize=false;
	public boolean printQuerySize=false;
	public boolean printContamHits=false;
	
	//For format 3
	public boolean useTaxidName=false;
	public boolean useImgName=false;
	public boolean useTaxName=false;
	
	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/
	
	public static final int default_maxRecords=20;
	public static final float default_minWKID=0.0001f;
	public static final int default_format=2;
	public static final boolean default_printTax=false;
	public static final boolean default_printOriginalName=false;
	public static final boolean default_printFileName=false;
	public static final boolean default_printImg=false;
	public static final boolean default_printAni=true;
	public static final boolean default_printCompleteness=true;
	public static final boolean default_printScore=false;
	public static final boolean default_printContam=true;
	
	public static final boolean default_printMatches=true;
	public static final boolean default_printLength=false;
	public static final boolean default_printTaxID=true;
	public static final boolean default_printGSize=true;
	public static final boolean default_printGKmers=false;
	public static final boolean default_printTaxName=true;
	public static final boolean default_printGSeqs=true;
	public static final boolean default_printGBases=false;

	public static final boolean default_printUMatches=true;
	public static final boolean default_printUContam=false;
	public static final boolean default_printNoHit=true;

	public static final boolean default_printColors=true;
	public static final int default_colorLevel=TaxTree.FAMILY_E;
	
	public static final int default_minHits=3;
	public static final int default_taxLevel=TaxTree.SPECIES;
	public static final int default_mode=SketchObject.ONE_SKETCH;
	public static final float default_samplerate=1;
	public static final long default_reads=-1;
	public static final int default_minKeyOccuranceCount=1;
	
	
	
}
