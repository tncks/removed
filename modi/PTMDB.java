package modi;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.StringTokenizer;
import java.util.TreeMap;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

import msutil.MSMass;

class GapKey extends Object {
	private Sequence seq;
	private int massDiff;
	private static int currentID = 0;
	private static final int unitSize = 10;
	private PTMPosition pos;
	
	public GapKey( Sequence seq, double massDiff, PTMPosition pos )
	{
		this.seq = seq;
		this.massDiff = (int)(massDiff*unitSize); // should be considered later. 4 just means a big enough tolerance for 0.5
		this.pos = pos;
	}
	
	public String	getString()		{ return seq.toString(); }
	public Sequence	getSeq()		{ return seq; }
	public double	getMassDiff()	{ return massDiff/(double)unitSize; }
	public int hashCode()
	{
		return seq.hashCode() + massDiff + pos.hashCode() ;
	}
	
	public boolean equals(Object key)
	{
		if (key instanceof GapKey)
			return this.seq.equals( ((GapKey)key).seq ) && this.massDiff==((GapKey)key).massDiff 
				&& this.pos==((GapKey)key).pos;
		else return false;
	}
}


public class PTMDB extends ArrayList<PTM> {
	private ArrayList<PTM>[][]		PTMTable;
	private TreeMap<Integer, String> classifications = new TreeMap<Integer, String>();
	public int parseUnimod() throws IOException
	{
		return parseUnimod(new File(Constants.UNIMOD_FILE_NAME));
	}
	public int parseUnimod(File unimodFile) throws IOException
	{
		System.out.print("[PTMDB]Parsing Unimod.xml");

		// Make document(DOM) from unimod.xml
        Document document; 
        DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
        try {
           DocumentBuilder builder = factory.newDocumentBuilder();
           document = builder.parse( unimodFile ); 
        } catch (SAXException sxe) {
           // Error generated during parsing)
           Exception  x = sxe;
           if (sxe.getException() != null)
               x = sxe.getException();
           x.printStackTrace();
           return -1;
        } catch (ParserConfigurationException pce) {
            // Parser with specified options can't be built
            pce.printStackTrace();
            return -1;
        } catch (IOException ioe) {
           // I/O error
           ioe.printStackTrace();
           return -1;
        }
       	
        // parse modifications
		Element modifications = (Element)document.getElementsByTagName("modifications").item(0);
		NodeList list = modifications.getElementsByTagName("modifications_row");
		class Modification {
			String		fullName;
			String		codeName;
			double		monoMass;
			double		avgeMass;
		}
		
		TreeMap<Integer, Modification>	modMap = new TreeMap<Integer, Modification>();
		for (int i=0; i<list.getLength(); i++)
		{
			Modification mod = new Modification();
			Element e = (Element)list.item(i);
			
			mod.fullName	= e.getAttribute("full_name");
			
			if( (mod.codeName= e.getAttribute("ex_code_name")).length() == 0 )
				mod.codeName	= e.getAttribute("code_name");
			mod.monoMass	= Double.parseDouble( e.getAttribute("mono_mass") );
			mod.avgeMass	= Double.parseDouble( e.getAttribute("avge_mass") );
			
			int modKey = Integer.parseInt(e.getAttribute("record_id"));
		
			modMap.put(modKey, mod);
		}
		
        // parse classifications
        Element classElement = (Element)document.getElementsByTagName("classifications").item(0);
        list = classElement.getElementsByTagName("classifications_row");
       	for(int i=0; i<list.getLength(); i++)
       	{
       		Element e = (Element)list.item(i);
       		classifications.put(Integer.parseInt(e.getAttribute("record_id")), e.getAttribute("classification"));
       	}
	
		// parse specificity & make PTMs
		PTMPosition[] arrPTMPosition = 
		{
			PTMPosition.NOT_DEFINED,
			PTMPosition.ANYWHERE,
			PTMPosition.ANY_N_TERM,
			PTMPosition.ANY_C_TERM,
			PTMPosition.PROTEIN_N_TERM,
			PTMPosition.PROTEIN_C_TERM
		};
		
		int id = 0;
		list = ((Element)document.getElementsByTagName("specificity").item(0)).getElementsByTagName("specificity_row");
        for (int i=0; i<list.getLength(); i++)
        {
        	Element 		e 		= (Element)list.item(i);
        	int 			modKey	= Integer.parseInt(e.getAttribute("mod_key"));
        	Modification	mod		= modMap.get(modKey);
        	String			residue	= e.getAttribute("one_letter");
        	boolean			residueNTerm = (residue.compareTo("N-term")==0);
        	boolean			residueCTerm = (residue.compareTo("C-term")==0);
        	AminoAcid		aa		= null;
        	PTMPosition		position = arrPTMPosition[Integer.parseInt(e.getAttribute("position_key"))-1];
        	String			classification = classifications.get(Integer.parseInt(e.getAttribute("classifications_key")));
        	
        	if (!residueNTerm && !residueCTerm) aa = AminoAcid.getAminoAcid(residue.charAt(0));
        	
    		// for unacceptable PTMs in Unimod.xml 
        	if (residueNTerm && position==PTMPosition.ANYWHERE) position = PTMPosition.ANY_N_TERM;
        	if (residueCTerm && position==PTMPosition.ANYWHERE) position = PTMPosition.ANY_C_TERM;
        	
        	
        	PTM ptm = new PTM(
        				id++,
        				mod.codeName,
        				mod.fullName,
        				mod.monoMass,
        				mod.avgeMass,
        				aa,
        				position,
        				classification
        			);
        	
        	this.add(ptm);
        }
        
		System.out.println( " -- [PTM DB] Entries : " + this.size() );
		return constructPTMTable();
	}

	public int parseUnimodOld() throws IOException	// deprecated, for parsing old version 
	{
		System.out.print("[PTMDB]Parsing Unimod.xml & exptmList.txt..");
		HashSet<String> exptmSet = new HashSet<String>();

		// Parse extpmList.txt
		BufferedReader in;
		try {
			in = new BufferedReader( new FileReader("exptmList.txt") );
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			return -1;
		}
		
		while (true)
		{
			String s = in.readLine();
			if (s==null) break;
			exptmSet.add(s);
		}
		in.close();		
		
		
		// Make document(DOM) from unimod.xml
        Document document; 
        DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
        try {
           DocumentBuilder builder = factory.newDocumentBuilder();
           document = builder.parse( new File(Constants.UNIMOD_FILE_NAME) ); 
        } catch (SAXException sxe) {
           // Error generated during parsing)
           Exception  x = sxe;
           if (sxe.getException() != null)
               x = sxe.getException();
           x.printStackTrace();
           return -1;
        } catch (ParserConfigurationException pce) {
            // Parser with specified options can't be built
            pce.printStackTrace();
            return -1;
        } catch (IOException ioe) {
           // I/O error
           ioe.printStackTrace();
           return -1;
        }
        
        // parse modifications
		Element modifications = (Element)document.getElementsByTagName("modifications").item(0);
		NodeList list = modifications.getElementsByTagName("row");
		class Modification {
			String		fullName;
			String		codeName;
			double		monoMass;
			double		avgeMass;
		}
		
		TreeMap<Integer, Modification>	modMap = new TreeMap<Integer, Modification>();;
		for (int i=0; i<list.getLength(); i++)
		{
			Modification mod = new Modification();
			Element e = (Element)list.item(i);
			
			mod.fullName	= e.getElementsByTagName("full_name").item(0).getTextContent();
			mod.codeName	= e.getElementsByTagName("code_name").item(0).getTextContent();
			mod.monoMass	= Double.parseDouble( e.getElementsByTagName("mono_mass").item(0).getTextContent() );
			mod.avgeMass	= Double.parseDouble( e.getElementsByTagName("avge_mass").item(0).getTextContent() );
			
			int modKey = Integer.parseInt(e.getElementsByTagName("record_id").item(0).getTextContent());
		
			modMap.put(modKey, mod);
		}
		
		// parse specificity & make PTMs
		PTMPosition[] arrPTMPosition = 
		{
			PTMPosition.NOT_DEFINED,
			PTMPosition.ANYWHERE,
			PTMPosition.ANY_N_TERM,
			PTMPosition.ANY_C_TERM,
			PTMPosition.PROTEIN_N_TERM,
			PTMPosition.PROTEIN_C_TERM
		};
		
		int id = 0;
		list = ((Element)document.getElementsByTagName("specificity").item(0)).getElementsByTagName("row");
        for (int i=0; i<list.getLength(); i++)
        {
        	Element 		e 		= (Element)list.item(i);
        	int 			modKey	= Integer.parseInt(e.getElementsByTagName("mod_key").item(0).getTextContent());
        	Modification	mod		= modMap.get(modKey);
        	String			residue	= e.getElementsByTagName("one_letter").item(0).getTextContent();
        	boolean			residueNTerm = (residue.compareTo("N-term")==0);
        	boolean			residueCTerm = (residue.compareTo("C-term")==0);
        	AminoAcid		aa		= null;
        	PTMPosition		position = arrPTMPosition[Integer.parseInt(e.getElementsByTagName("position_key").item(0).getTextContent())-1]; 
        	
        	if (exptmSet.contains(mod.codeName)) 
        		continue; // Skip ptm in exptList.txt
        	
        	if (!residueNTerm && !residueCTerm) aa = AminoAcid.getAminoAcid(residue.charAt(0));
        	
    		// for unacceptable PTMs in Unimod.xml 
        	if (residueNTerm && position==PTMPosition.ANYWHERE) position = PTMPosition.ANY_N_TERM;
        	if (residueCTerm && position==PTMPosition.ANYWHERE) position = PTMPosition.ANY_C_TERM;
        	
        	
        	PTM ptm = new PTM(
        				id++,
        				mod.codeName,
        				mod.fullName,
        				mod.monoMass,
        				mod.avgeMass,
        				aa,
        				position
        			);
        	
        	this.add(ptm);
        }
        
		System.out.println("Done!! Entries: "+this.size());
		return constructPTMTable();
	}
	
	public int parseUserDB(String fileName) throws Exception
	{
		System.out.print("[PTMDB]Parsing User DB..");
		org.jdom.Document document = null;
		try {
			 document = new org.jdom.input.SAXBuilder().build(new File(fileName));
		} catch(org.jdom.JDOMException e) {
			System.out.println(e);
		}
		
		org.jdom.Element rootElement = document.getRootElement();
		org.jdom.Element classificationsElement = rootElement.getChild("classifications");
		if(classificationsElement != null)	// for backward compatibility
		{
			classifications.clear();
			for(Object obj : classificationsElement.getChildren())
			{
				org.jdom.Element classificationRowElement = (org.jdom.Element)obj;
				classifications.put(Integer.parseInt(classificationRowElement.getAttributeValue("recordID")),
						classificationRowElement.getAttributeValue("classification"));
			}
		}		
		int id=0;
		for(Object obj : rootElement.getChildren("PTM"))
		{
			org.jdom.Element elemPTM = (org.jdom.Element)obj;

			String residueStr = elemPTM.getChildText("residue");
			AminoAcid residue = null;
			if(residueStr != null && residueStr.compareToIgnoreCase("N-term") != 0 && residueStr.compareToIgnoreCase("C-term") != 0)
			{
				assert(residueStr.length() == 1);
				residue = AminoAcid.getAminoAcid(residueStr.charAt(0));
			}
			
			PTMPosition position = null;
			String positionStr = elemPTM.getChildText("position");
			if(positionStr.equalsIgnoreCase("ANYWHERE"))
				position = PTMPosition.ANYWHERE;
			else if(positionStr.equalsIgnoreCase("ANY_N_TERM"))
				position = PTMPosition.ANY_N_TERM;
			else if(positionStr.equalsIgnoreCase("ANY_C_TERM"))
				position = PTMPosition.ANY_C_TERM;
			else if(positionStr.equalsIgnoreCase("PROTEIN_N_TERM"))
				position = PTMPosition.PROTEIN_N_TERM;
			else if(positionStr.equalsIgnoreCase("PROTEIN_C_TERM"))
				position = PTMPosition.PROTEIN_C_TERM;
			else
				assert(false);

			double shift = ( residue != null && residue.getIndex() == 1 )? Constants.alkylatedToCys : 0;
			this.add(new PTM(
						id++,
						elemPTM.getChildText("name"), 
						elemPTM.getChildText("fullName"), 
						Double.parseDouble(elemPTM.getChildText("massDifference"))-shift,
						Double.parseDouble(elemPTM.getChildText("avgMassDifference"))-shift,
						residue,
						position
						)
					);
		}
		
		System.out.println("Done!! Entries: "+this.size());
		return constructPTMTable();
	}

	public boolean isNovelPtm(int site, int res, double mass){
		int i;
		for(i=0; i<PTMTable[res][PTMPosition.ANYWHERE.ordinal()].size(); i++ ){
			int delta= (int)Math.round(PTMTable[res][PTMPosition.ANYWHERE.ordinal()].get(i).getMassDifference());
			if( delta == (int)Math.round(mass) )
				return false;		
		}		
		if( site == 0 ){
			for(i=0; i<PTMTable[res][PTMPosition.ANY_N_TERM.ordinal()].size(); i++ ){
				int delta= (int)Math.round(PTMTable[res][PTMPosition.ANY_N_TERM.ordinal()].get(i).getMassDifference());
				if( delta == (int)Math.round(mass) )
					return false;	
			}
		}
		return true;	
	}
	
	public void constructPTMLookupTable(){
		constructPTMTable();
	}

	private int constructPTMTable() {	
		this.sortByPTMPosition();
				
		PTMTable = new ArrayList[AminoAcid.getIndexSize()][PTMPosition.PTMPOSITION_COUNT.ordinal()];
		for (int i=0; i<AminoAcid.getIndexSize(); i++) {
			for (int j=0; j<PTMPosition.PTMPOSITION_COUNT.ordinal(); j++)
				PTMTable[i][j] = new ArrayList<PTM>();		
		}
		
		
		for ( PTM ptm : this ) {
			if( ptm.getMassDifference() < Constants.minModifiedMass || ptm.getMassDifference() > Constants.maxModifiedMass ) continue;
			if ( ptm.getResidue() != null ) {
				int aa = ptm.getResidue().getIndex();
				PTMTable[ptm.getResidue().getIndex()][ptm.getPTMPosition().ordinal()].add(ptm);
				
				if( ptm.getPTMPosition() == PTMPosition.ANYWHERE ){
					if( modRange[0][aa] > ptm.getMassDifference() ) modRange[0][aa] = ptm.getMassDifference(); // Min
					else if( modRange[1][aa] < ptm.getMassDifference() ) modRange[1][aa] = ptm.getMassDifference(); // Max
				}
				else if( ptm.getPTMPosition() == PTMPosition.ANY_N_TERM || ptm.getPTMPosition() == PTMPosition.PROTEIN_N_TERM ){
					if( ntermModRange[0][aa] > ptm.getMassDifference() ) ntermModRange[0][aa] = ptm.getMassDifference(); // Min
					else if( ntermModRange[1][aa] < ptm.getMassDifference() ) ntermModRange[1][aa] = ptm.getMassDifference(); // Max
				}
				else if( ptm.getPTMPosition() == PTMPosition.ANY_C_TERM || ptm.getPTMPosition() == PTMPosition.PROTEIN_C_TERM ){
					if( ctermModRange[0][aa] > ptm.getMassDifference() ) ctermModRange[0][aa] = ptm.getMassDifference(); // Min
					else if( ctermModRange[1][aa] < ptm.getMassDifference() ) ctermModRange[1][aa] = ptm.getMassDifference(); // Max
				}
			}
			else {
				if ( !ptm.isResidueCTerm() && !ptm.isResidueNTerm() ) return -1;
				
				for (int j=0; j<AminoAcid.getIndexSize(); j++){				
					
					if( PTMTable[j][PTMPosition.ANYWHERE.ordinal()].contains(ptm) ) continue;
					
					if ( ptm.getPTMPosition() == PTMPosition.PROTEIN_N_TERM ) {
						if( PTMTable[j][PTMPosition.ANY_N_TERM.ordinal()].contains(ptm) )
							continue;
					}					
					if ( ptm.getPTMPosition() == PTMPosition.PROTEIN_C_TERM ) {
						if( PTMTable[j][PTMPosition.ANY_C_TERM.ordinal()].contains(ptm) )
							continue;
					}//*/

					PTMTable[j][ptm.getPTMPosition().ordinal()].add(ptm);
					
					if( ptm.getPTMPosition() == PTMPosition.ANY_N_TERM || ptm.getPTMPosition() == PTMPosition.PROTEIN_N_TERM ){
						if( ntermModRange[0][j] > ptm.getMassDifference() ) ntermModRange[0][j] = ptm.getMassDifference(); // Min
						else if( ntermModRange[1][j] < ptm.getMassDifference() ) ntermModRange[1][j] = ptm.getMassDifference(); // Max
					}
					else if( ptm.getPTMPosition() == PTMPosition.ANY_C_TERM || ptm.getPTMPosition() == PTMPosition.PROTEIN_C_TERM ){
						if( ctermModRange[0][j] > ptm.getMassDifference() ) ctermModRange[0][j] = ptm.getMassDifference(); // Min
						else if( ctermModRange[1][j] < ptm.getMassDifference() ) ctermModRange[1][j] = ptm.getMassDifference(); // Max
					}
				}
			}
		}
		
		Collections.sort(this);	
		return 0;
	}
	
	// Temporary valuables for recursive search
	private Sequence			seq;
	private	double				massDiff;
	private	PTMPosition			position;
	private ArrayList<PTMRun>	result;
	private PTM[]				occur;
	private int[]				numNextFixSite;
	private int					numMaxMods;

	// Finding PTMs by Exhaustive Depth-first-search
	private void findPTM_DFS( double mass, int pos, int cnt, int extra )
	{	
		// when search has reached the leaf node(end of sequence), check result
		if ( pos == seq.size() ) {			
			double ierror = Math.abs( mass - massDiff );
		//	System.out.println(Constants.gapTolerance);
		//	if ( ierror <= Constants.gapTolerance ) 
			if ( ierror <= Constants.gapTolerance && Constants.isWithinAccuracy(ierror) ) 
			{
				PTMRun run = new PTMRun();
				for (int i=0; i<seq.size(); i++)
					if ( occur[i] != null ){
						run.add(new PTMOccurrence(i, occur[i]));
					}
				
				if ( run.size() > 0 ) {
					run.setError( ierror );
					this.result.add(run);
				}
			}
			return;
		}
		
		// recursion call without checking PTM on this pos
		findPTM_DFS( mass, pos+1, cnt, extra );		
		
		if( cnt >= Constants.maxPTMPerGap ) return;
		if( extra == 1 && numMaxMods == 1 ) {
			if ( numNextFixSite[pos] == 0 ) return;
		}

		int residueIndex = seq.get(pos).getIndex();
		for (int i=1; i<PTMPosition.PTMPOSITION_COUNT.ordinal(); i++) {
			if (pos>0 && pos<seq.size()-1 && i!=1) continue;
			if ((i==2 || i==4) && pos!=0) continue;
			if ((i==3 || i==5) && pos!=seq.size()-1) continue;
			
			if ((i==2) && this.position!=PTMPosition.ANY_N_TERM 
					&& this.position!=PTMPosition.PROTEIN_N_TERM ) continue;
			if ((i==3) && this.position!=PTMPosition.ANY_C_TERM 
					&& this.position!=PTMPosition.PROTEIN_C_TERM ) continue;
			
			if ((i==4) && this.position!=PTMPosition.PROTEIN_N_TERM ) continue;
			if ((i==5) && this.position!=PTMPosition.PROTEIN_C_TERM ) continue;
			
			for (PTM ptm : PTMTable[residueIndex][i]) {				
				occur[pos] = ptm;				
				if( extra + ptm.getModCount() > numMaxMods ) continue;
				findPTM_DFS( mass + ptm.getMassDifference(), pos+1, cnt+1, extra+ptm.getModCount() );
						
			}
		}
		occur[pos] = null;
	}

	// Hash table for reusing the result for Gaps already searched once
	private HashMap< GapKey, PTMSearchResult >	hashTable = new HashMap< GapKey, PTMSearchResult >();
	
	public void clearHashTable()		{ hashTable.clear(); }
	
	// before ver 0.8, deprecated, Not Cashing
	public PTMSearchResult searchPTM( Sequence seq, double massDiff, PTMPosition position ) {		
		PTMSearchResult searchResult;
		// Hash miss
		ArrayList<PTMRun> newGapInterpret = new ArrayList<PTMRun>();
	
		double ierror = Math.abs( massDiff );
		if( ierror < Constants.nonModifiedDelta ){
			searchResult = new PTMSearchResult( newGapInterpret, true );			
			return searchResult;
		}
		
		if( ierror < Constants.gapTolerance ){
			PTMRun run = new PTMRun();
			run.setError(ierror);
			newGapInterpret.add( run );			
		}
		
		// PTM search
		this.seq 		= seq;
		this.massDiff	= massDiff;
		this.position	= position;
		this.result		= newGapInterpret;
		occur = new PTM[seq.size()];
		numNextFixSite = new int[seq.size()];				
		numMaxMods = Constants.getMaxPTMOccurrence(seq.size());
				
		int cum = ( Constants.CTERM_FIX_MOD !=0 && (this.position == PTMPosition.ANY_C_TERM || this.position == PTMPosition.PROTEIN_C_TERM) )? 1 : 0;
		for(int i=seq.size()-1; i>-1; i--){
			if( seq.get(i).isLabelled() ) {
				numNextFixSite[i] = ++cum;
			}
			else {				
				numNextFixSite[i] = cum;
			}
		}
		
		findPTM_DFS( 0.0, 0, 0, 0 );
		this.result		= null;		// release reference from PTMDB class
		
		searchResult = new PTMSearchResult(newGapInterpret, newGapInterpret != null && newGapInterpret.size() > 0);
		
	//	hashTable.put( new GapKey(seq, massDiff, position), searchResult );
		
		return searchResult;
	}
	
	public org.jdom.Element getGapInterpretListElement()
	{
		if(Constants.ANALYSIS_VERSION < 1.0)
			return null;
		
		org.jdom.Element gapInterpretListElement = new org.jdom.Element("gapInterpretList");
		// HashMap< GapKey, ArrayList<PTMRun> > hashTable
		Iterator<Map.Entry<GapKey, PTMSearchResult>> it = hashTable.entrySet().iterator();
		
		// to print actualInterpret entry by order of id
		
		HashMap<Integer, GapKey> gapKeyOrder = new HashMap<Integer, GapKey>();
		while(it.hasNext())
		{
			Map.Entry<GapKey, PTMSearchResult> entry = it.next();
			gapKeyOrder.put(entry.getValue().getID(), entry.getKey());
		}

		ArrayList<Integer> orderedIDSet = new ArrayList<Integer> (gapKeyOrder.keySet());
		Collections.sort(orderedIDSet);
			
		for(Integer keyID : orderedIDSet)
		{
			GapKey key = gapKeyOrder.get(keyID);
			PTMSearchResult result = hashTable.get(key);
			
			org.jdom.Element elemGapInterpret = new org.jdom.Element("gapInterpret");
			elemGapInterpret.setAttribute(new org.jdom.Attribute("id", String.valueOf(result.getID())));
			elemGapInterpret.setAttribute(new org.jdom.Attribute("sequence", key.getString()));
			elemGapInterpret.setAttribute(new org.jdom.Attribute("offset", String.valueOf(key.getMassDiff())));
			String hasPTMStr;
			ArrayList<PTMRun> ptmRunList = result.getPTMRun();
			if(ptmRunList == null || ptmRunList.isEmpty())
				hasPTMStr = "no";
			else
				hasPTMStr = "yes";
			elemGapInterpret.setAttribute(new org.jdom.Attribute("hasPTM", hasPTMStr));
			int id = 0;
			for(PTMRun ptmRun : ptmRunList)
			{
				org.jdom.Element elemActualInterpret = new org.jdom.Element("actualInterpret");
				elemActualInterpret.setAttribute(new org.jdom.Attribute("id", String.valueOf(id++)));
				
				org.jdom.Element elemPTMRun = new org.jdom.Element("PTMRun");
				for(PTMOccurrence occr : ptmRun)
				{
					org.jdom.Element elemPTMOccr = new org.jdom.Element("PTMOccr");
					elemPTMOccr.setAttribute(new org.jdom.Attribute("id", String.valueOf(occr.getPTM().getID())));
					elemPTMOccr.setAttribute(new org.jdom.Attribute("position", String.valueOf(occr.getPosition())));
					elemPTMRun.addContent(elemPTMOccr);
				}
				elemActualInterpret.addContent(elemPTMRun);
				
				org.jdom.Element elemBTheoPeaks = new org.jdom.Element("bTheoPeaks");
				elemBTheoPeaks.setText(Gap.getTheoreticalMassStr(getBTheoreticalPeaks(key.getSeq(), ptmRun)));
				elemActualInterpret.addContent(elemBTheoPeaks);

				org.jdom.Element elemYTheoPeaks = new org.jdom.Element("yTheoPeaks");
				elemYTheoPeaks.setText(Gap.getTheoreticalMassStr(getYTheoreticalPeaks(key.getSeq(), ptmRun)));
				elemActualInterpret.addContent(elemYTheoPeaks);
				
				elemGapInterpret.addContent(elemActualInterpret);
			}
			gapInterpretListElement.addContent(elemGapInterpret);
		}
		return gapInterpretListElement;
	}

	public	ArrayList<Peak> getBTheoreticalPeaks(Sequence seq, PTMRun run)
	{
		double [] ptmMass = new double[seq.size()];
		for(PTMOccurrence occr : run)
			ptmMass[occr.getPosition()] += occr.getPTM().getMassDifference();

		ArrayList<Peak> bPeaks = new ArrayList<Peak>();
		double mass = 0;
		for(int i=0; i<seq.size()-1; i++)
		{
			mass = mass + seq.get(i).getMass() + ptmMass[i];
			bPeaks.add(new Peak(-1, mass, 0.));
		}
		return bPeaks;
	}
	
	public	ArrayList<Peak> getYTheoreticalPeaks(Sequence seq, PTMRun run)
	{
		double [] ptmMass = new double[seq.size()];
		for(PTMOccurrence occr : run)
			ptmMass[occr.getPosition()] += occr.getPTM().getMassDifference();

		ArrayList<Peak> yPeaks = new ArrayList<Peak>();
		double mass = 0;
		for(int i=0; i<seq.size()-1; i++)
		{
			mass = mass + seq.get(seq.size()-1-i).getMass() + ptmMass[i];
			yPeaks.add(new Peak(-1, mass, 0.));
		}
		return yPeaks;
	}
	
	public	org.jdom.Element getElement()
	{
		org.jdom.Element elemPTMDB = new org.jdom.Element("PTMDB");
		for(PTM ptm : this)
		{
			org.jdom.Element elemPTM = new org.jdom.Element("PTM");
			elemPTM.setAttribute(new org.jdom.Attribute("id", String.valueOf(ptm.getID())));
			elemPTM.setAttribute(new org.jdom.Attribute("name", ptm.getName()));
			elemPTM.setAttribute(new org.jdom.Attribute("fullName", ptm.getFullName()));
			elemPTM.setAttribute(new org.jdom.Attribute("massDiff", Constants.getString(ptm.getMassDifference())));
			if(ptm.getResidue() == null)
				elemPTM.setAttribute(new org.jdom.Attribute("residue", "ANY"));
			else
				elemPTM.setAttribute(new org.jdom.Attribute("residue", ptm.getResidue().toString()));
			elemPTM.setAttribute(new org.jdom.Attribute("position", ptm.getPTMPosition().toString()));
			elemPTMDB.addContent(elemPTM);
		}
		
		return elemPTMDB;
	}
	// for PTMDB.xml output
	public	org.jdom.Document getDocument()
	{
		org.jdom.Document document = new org.jdom.Document();
		org.jdom.Element elemPTMDB = new org.jdom.Element("PTMDB");
		elemPTMDB.setAttribute(new org.jdom.Attribute("majorVersion", "1"));
		elemPTMDB.setAttribute(new org.jdom.Attribute("minorVersion", "0"));
		document.addContent(elemPTMDB);

		org.jdom.Element elemClassifications = new org.jdom.Element("classifications");
		elemPTMDB.addContent(elemClassifications);

		Iterator<Map.Entry<Integer, String>> it = classifications.entrySet().iterator();
		while(it.hasNext())
		{
			org.jdom.Element elemClassificationsRow = new org.jdom.Element("classificationRow");
			Map.Entry<Integer, String> entry = it.next();
			elemClassificationsRow.setAttribute(new org.jdom.Attribute("recordID", entry.getKey().toString()));
			elemClassificationsRow.setAttribute(new org.jdom.Attribute("classification", entry.getValue()));
			elemClassifications.addContent(elemClassificationsRow);
		}
		
		for(PTM ptm : this)
		{
			
		//	if( ptm.getClassification().contains("label") ) continue;
			
			org.jdom.Element elemPTM = new org.jdom.Element("PTM");
			elemPTMDB.addContent(elemPTM);
			
			org.jdom.Element elemName = new org.jdom.Element("name");
			elemName.setText(ptm.getName());
			elemPTM.addContent(elemName);

			org.jdom.Element elemFullName = new org.jdom.Element("fullName");
			elemFullName.setText(ptm.getFullName());
			elemPTM.addContent(elemFullName);

			org.jdom.Element elemClassification = new org.jdom.Element("classification");
			elemClassification.setText(ptm.getClassification());
			elemPTM.addContent(elemClassification);
			
			org.jdom.Element elemMassDifference = new org.jdom.Element("massDifference");
			elemMassDifference.setText(String.valueOf(ptm.getMassDifference()));
			elemPTM.addContent(elemMassDifference);

			org.jdom.Element elemAvgMassDifference = new org.jdom.Element("avgMassDifference");
			elemAvgMassDifference.setText(String.valueOf(ptm.getAvgMassDifference()));
			elemPTM.addContent(elemAvgMassDifference);
			
			org.jdom.Element elemResidue = new org.jdom.Element("residue");
			if(ptm.getResidue() != null)
				elemResidue.setText(ptm.getResidue().toString());
			else
			{
				if(ptm.isResidueNTerm()) 
					elemResidue.setText("N-term");
				else if(ptm.isResidueCTerm())
					elemResidue.setText("C-term");
				else
					assert(false);
			}
			
				
			elemPTM.addContent(elemResidue);

			org.jdom.Element elemPosition = new org.jdom.Element("position");
			elemPosition.setText(ptm.getPTMPosition().toString());
			elemPTM.addContent(elemPosition);
		}
		
		return document;
	}

	// Since ver 1.0
	public	org.jdom.Document xxxxxxgetDocumentNew()
	{
		org.jdom.Document document = new org.jdom.Document();
		org.jdom.Element elemPTMDB = new org.jdom.Element("PTMDB");
		elemPTMDB.setAttribute(new org.jdom.Attribute("majorVersion", "1"));
		elemPTMDB.setAttribute(new org.jdom.Attribute("minorVersion", "0"));
		document.addContent(elemPTMDB);

		org.jdom.Element elemClassifications = new org.jdom.Element("classifications");
		elemPTMDB.addContent(elemClassifications);

		Iterator<Map.Entry<Integer, String>> it = classifications.entrySet().iterator();
		while(it.hasNext())
		{
			org.jdom.Element elemClassificationsRow = new org.jdom.Element("classificationRow");
			Map.Entry<Integer, String> entry = it.next();
			elemClassificationsRow.setAttribute(new org.jdom.Attribute("recordID", entry.getKey().toString()));
			elemClassificationsRow.setAttribute(new org.jdom.Attribute("classification", entry.getValue()));
			elemClassifications.addContent(elemClassificationsRow);
		}
		
		for(PTM ptm : this)
		{
			org.jdom.Element elemPTM = new org.jdom.Element("PTM");
			elemPTMDB.addContent(elemPTM);
			
			org.jdom.Element elemName = new org.jdom.Element("name");
			elemName.setText(ptm.getName());
			elemPTM.addContent(elemName);

			org.jdom.Element elemFullName = new org.jdom.Element("fullName");
			elemFullName.setText(ptm.getFullName());
			elemPTM.addContent(elemFullName);

			org.jdom.Element elemClassification = new org.jdom.Element("classification");
			elemClassification.setText(ptm.getClassification());
			elemPTM.addContent(elemClassification);

			org.jdom.Element elemMassDifference = new org.jdom.Element("massDifference");
			elemMassDifference.setText(Constants.getString(ptm.getMassDifference()));
			elemPTM.addContent(elemMassDifference);

			org.jdom.Element elemAvgMassDifference = new org.jdom.Element("avgMassDifference");
			elemAvgMassDifference.setText(Constants.getString(ptm.getAvgMassDifference()));
			elemPTM.addContent(elemAvgMassDifference);
			
			org.jdom.Element elemResidue = new org.jdom.Element("residue");
			if(ptm.getResidue() != null)
				elemResidue.setText(ptm.getResidue().toString());
			else
			{
				if(ptm.isResidueNTerm()) 
					elemResidue.setText("N-term");
				else if(ptm.isResidueCTerm())
					elemResidue.setText("C-term");
				else
					assert(false);
			}
			
				
			elemPTM.addContent(elemResidue);

			org.jdom.Element elemPosition = new org.jdom.Element("position");
			elemPosition.setText(ptm.getPTMPosition().toString());
			elemPTM.addContent(elemPosition);
		}
		
		return document;
	}

	public	PTM	getPTM(String name, AminoAcid residue) {
		for(PTM ptm : this)
			if(ptm.getName().equalsIgnoreCase(name) && ptm.getResidue() == residue)
				return ptm;
		return null;
	}
	
	public PTM	getPTM(String name, char residue){	
		int AA = AminoAcid.getAminoAcid(residue).getIndex();
		for(int j=0; j<PTMPosition.PTMPOSITION_COUNT.ordinal(); j++){
			for ( PTM ptm : PTMTable[AA][j] ){
				if( ptm.getName().equals(name) ){
					return ptm;
				}
			}
		}
		return null;
	}
	
	public void sortByPTMPosition(){
		Collections.sort( this, new PTMPosComparator() );
	}

	public String toPRIX(){
		StringBuffer des = new StringBuffer(this.size()+"|");
		for( PTM p : this ){
			des.append( p.getDescription() );
		}
		return des.toString();
	}
	public String toPRIX(String insertMod, String AA, double diff){
		StringBuffer des = new StringBuffer( (this.size()+1)+"|" );
		des.append( "*|C|-1|" );
		for( PTM p : this ){
			des.append(p.getDescription());
		}
		return des.toString();
	}
	
	public void toPRIXParam(String outputName) throws IOException{
		PrintWriter ptms = new PrintWriter(new BufferedWriter(new FileWriter(outputName+".ptms")));
		for( PTM p : this ){
			ptms.printf( "<mod name=\"%s\" site=\"%s\" position=\"%s\" massdiff=\"%f\" />\r\n", 
					p.getName(), p.getSite(), p.getPTMPosition(), p.getMassDifference() );
		}
		ptms.close();
	}
	
	public int setVariableModificatinos( String fileName, double[] fixedAA, boolean canBeModifiedOnFixedAA ) throws Exception
	{
		int id = this.size();
		org.jdom.Document document = null;
		try {
			 document = new org.jdom.input.SAXBuilder().build(new File(fileName));
		} catch(org.jdom.JDOMException e) {
			System.out.println(e);
		}
		
		org.jdom.Element rootElement = document.getRootElement();
		org.jdom.Element classificationsElement = rootElement.getChild("classifications");
		if(classificationsElement != null)	// for backward compatibility
		{
			classifications.clear();
			for(Object obj : classificationsElement.getChildren())
			{
				org.jdom.Element classificationRowElement = (org.jdom.Element)obj;
				classifications.put(Integer.parseInt(classificationRowElement.getAttributeValue("recordID")),
						classificationRowElement.getAttributeValue("classification"));
			}
		}		

		for(Object obj : rootElement.getChildren("PTM"))
		{
			org.jdom.Element elemPTM = (org.jdom.Element)obj;

			String pname = elemPTM.getChildText("name");
			
			String category = elemPTM.getChildText("classification");
			
			String residueStr = elemPTM.getChildText("residue");
			AminoAcid residue = null;
			
			double massdelta = 0;
			if(residueStr != null && residueStr.compareToIgnoreCase("N-term") != 0 && residueStr.compareToIgnoreCase("C-term") != 0) {
				assert( residueStr.length() == 1 );
				residue = AminoAcid.getAminoAcid(residueStr.charAt(0));
				if( fixedAA[residueStr.charAt(0)-'A'] != 0 ) {
					if( canBeModifiedOnFixedAA ) massdelta = fixedAA[residueStr.charAt(0)-'A'];
					else continue;
				}
				if( "AA substitution".compareTo(category) == 0 ){
					char tarAA = AminoAcid.getAminoAcidBy3Letter(pname.substring(pname.length()-3));
					massdelta -= fixedAA[tarAA-'A'];
					if( tarAA == 'I' && fixedAA['I'-'A'] != fixedAA['L'-'A'] ){
						System.out.println("PROGRAM SHOULD BE FIXED");
						System.exit(1);
					}
				}
			}		
			
			PTMPosition position = null;
			String positionStr = elemPTM.getChildText("position");
			if(positionStr.equalsIgnoreCase("ANYWHERE"))
				position = PTMPosition.ANYWHERE;
			else if(positionStr.equalsIgnoreCase("ANY_N_TERM"))
				position = PTMPosition.ANY_N_TERM;
			else if(positionStr.equalsIgnoreCase("ANY_C_TERM"))
				position = PTMPosition.ANY_C_TERM;
			else if(positionStr.equalsIgnoreCase("PROTEIN_N_TERM"))
				position = PTMPosition.PROTEIN_N_TERM;
			else if(positionStr.equalsIgnoreCase("PROTEIN_C_TERM"))
				position = PTMPosition.PROTEIN_C_TERM;
			else
				assert(false);

		//	if( Constants.NTERM_FIX_MOD != 0 && (position == PTMPosition.ANY_N_TERM || position == PTMPosition.PROTEIN_N_TERM) ) continue;
		//	if( Constants.CTERM_FIX_MOD != 0 && (position == PTMPosition.ANY_C_TERM || position == PTMPosition.PROTEIN_C_TERM) ) continue;

			if( Constants.NTERM_FIX_MOD != 0 && (position == PTMPosition.ANY_N_TERM || position == PTMPosition.PROTEIN_N_TERM) ) {
				if( canBeModifiedOnFixedAA ) {
					pname += "/Nterm";
					massdelta += Constants.NTERM_FIX_MOD;
				}
				else continue;
			}
			if( Constants.CTERM_FIX_MOD != 0 && (position == PTMPosition.ANY_C_TERM || position == PTMPosition.PROTEIN_C_TERM) ) {
				if( canBeModifiedOnFixedAA ) {
					massdelta += Constants.CTERM_FIX_MOD;
					pname += "/Cterm";
				}
				else continue;
			}
			
			double ac_delta = Double.parseDouble(elemPTM.getChildText("massDifference"));
			if( ac_delta < Constants.minModifiedMass || ac_delta > Constants.maxModifiedMass ) continue;
			ac_delta -= massdelta;
			
			this.add(new PTM(
						id++,
						pname, 
						"", 
						ac_delta,
						0,
						residue,
						position,
						category
						)
					);
		}
		return 1;
	}
	
	public int setVariableModificatinos( org.jdom.Element modifications, double[] fixedAA, boolean canBeModifiedOnFixedAA ) throws Exception
	{
		int id = this.size();
		for(Object obj : modifications.getChildren("mod"))
		{
			org.jdom.Element elemPTM = (org.jdom.Element)obj;
			
			String residueStr = elemPTM.getAttributeValue("site");
			AminoAcid residue = null;
			double massdelta = 0;
			if(residueStr != null && residueStr.compareToIgnoreCase("N-term") != 0 && residueStr.compareToIgnoreCase("C-term") != 0) {
				assert(residueStr.length() == 1);
				residue = AminoAcid.getAminoAcid(residueStr.charAt(0));
				if( fixedAA[residueStr.charAt(0)-'A'] != 0 ) {
					if( canBeModifiedOnFixedAA ) massdelta = fixedAA[residueStr.charAt(0)-'A'];
					else continue;
				}
			}			
			
			PTMPosition position = null;
			String positionStr = elemPTM.getAttributeValue("position");
			if(positionStr.equalsIgnoreCase("ANYWHERE"))
				position = PTMPosition.ANYWHERE;
			else if(positionStr.equalsIgnoreCase("ANY_N_TERM"))
				position = PTMPosition.ANY_N_TERM;
			else if(positionStr.equalsIgnoreCase("ANY_C_TERM"))
				position = PTMPosition.ANY_C_TERM;
			else if(positionStr.equalsIgnoreCase("PROTEIN_N_TERM"))
				position = PTMPosition.PROTEIN_N_TERM;
			else if(positionStr.equalsIgnoreCase("PROTEIN_C_TERM"))
				position = PTMPosition.PROTEIN_C_TERM;
			else
				assert(false);

			String pname = elemPTM.getAttributeValue("name");
			
			if( Constants.NTERM_FIX_MOD != 0 && (position == PTMPosition.ANY_N_TERM || position == PTMPosition.PROTEIN_N_TERM) ) {
				if( canBeModifiedOnFixedAA ) {
					pname += "/Nterm";
					massdelta += Constants.NTERM_FIX_MOD;
				}
				else continue;
			}
			if( Constants.CTERM_FIX_MOD != 0 && (position == PTMPosition.ANY_C_TERM || position == PTMPosition.PROTEIN_C_TERM) ) {
				if( canBeModifiedOnFixedAA ) {
					pname += "/Cterm";
					massdelta += Constants.CTERM_FIX_MOD;
				}
				else continue;
			}
			
			double ac_delta = Double.parseDouble(elemPTM.getAttributeValue("massdiff"));
			if( ac_delta < Constants.minModifiedMass || ac_delta > Constants.maxModifiedMass ) continue;
			ac_delta -= massdelta;
			
			this.add(new PTM(
						id++,
						pname, 
						"", 
						ac_delta,
						0,
						residue,
						position
						)
					);
		}
		return 1;//constructPTMTable();
	}
	
	public int setVariableModificatinos( ArrayList<String> bondPtms, double[] fixedMod ){//DBOND
		int id = this.size();					
		for( String s : bondPtms ){
			//name  where  residue  mass
			StringTokenizer token = new StringTokenizer(s, ", \t");
			if( token.countTokens() == 4 ){
				String name = token.nextToken();
				String where = token.nextToken();
				String res = token.nextToken();
				double mass = Double.parseDouble(token.nextToken());

				PTMPosition position = null;
				if( where.toUpperCase().equals("NTERM") )
					position = PTMPosition.ANY_N_TERM;
				else if( where.toUpperCase().equals("CTERM") )
					position = PTMPosition.ANY_C_TERM;
				else
					position = PTMPosition.ANYWHERE;
				
				AminoAcid residue = null;
				double massdelta = 0;
				if( !res.toUpperCase().equals("NTERM") && !res.toUpperCase().equals("CTERM") ){
					residue = AminoAcid.getAminoAcid(res.charAt(0));
					if( fixedMod[res.charAt(0)-'A'] != 0 )  massdelta = fixedMod[res.charAt(0)-'A'];
				}
				
				if( Constants.NTERM_FIX_MOD != 0 && (position == PTMPosition.ANY_N_TERM || position == PTMPosition.PROTEIN_N_TERM) ) {
					massdelta += Constants.NTERM_FIX_MOD;
				}
				if( Constants.CTERM_FIX_MOD != 0 && (position == PTMPosition.ANY_C_TERM || position == PTMPosition.PROTEIN_C_TERM) ) {
					massdelta += Constants.CTERM_FIX_MOD;
				}
				
				this.add( new PTM(
								id++,
								name, 
								"", 
								mass-massdelta,
								0,
								residue,
								position
								) );			
			}
		}
		return constructPTMTable();
		
	}
	
	public void write( String fileName ) throws Exception
	{
		PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(fileName)));
		for( PTM p : this ){
			out.println(p.getName()+"\t"+p.getSite()+"\t"+p.getPTMPosition()+"\t"+p.getMassDifference());
		}
		out.close();
	}
	
	public int setFixedModificatinos(org.jdom.Element modifications, double[] fixedAA) throws Exception
	{	
		int id=0;
		for(Object obj : modifications.getChildren("mod"))
		{
			org.jdom.Element elemPTM = (org.jdom.Element)obj;
			
			String residueStr = elemPTM.getAttributeValue("site");
			AminoAcid residue = null;
			
			if( residueStr == null ) continue;
			else if( residueStr.compareToIgnoreCase("N-term") == 0 ){
				Constants.NTERM_FIX_MOD += Double.parseDouble(elemPTM.getAttributeValue("massdiff"));
			}
			else if( residueStr.compareToIgnoreCase("C-term") == 0 ){
				Constants.CTERM_FIX_MOD += Double.parseDouble(elemPTM.getAttributeValue("massdiff"));
			}
			else residue = AminoAcid.getAminoAcid(residueStr.charAt(0));
			
			double diff = Double.parseDouble(elemPTM.getAttributeValue("massdiff"));		
			
			PTMPosition position = null;
			String positionStr = elemPTM.getAttributeValue("position");
			if(positionStr.equalsIgnoreCase("ANYWHERE"))
				position = PTMPosition.ANYWHERE;
			else if(positionStr.equalsIgnoreCase("ANY_N_TERM"))
				position = PTMPosition.ANY_N_TERM;
			else if(positionStr.equalsIgnoreCase("ANY_C_TERM"))
				position = PTMPosition.ANY_C_TERM;
			else if(positionStr.equalsIgnoreCase("PROTEIN_N_TERM"))
				position = PTMPosition.PROTEIN_N_TERM;
			else if(positionStr.equalsIgnoreCase("PROTEIN_C_TERM"))
				position = PTMPosition.PROTEIN_C_TERM;
			else
				assert(false);

			this.add(new PTM(
						id++,
						elemPTM.getAttributeValue("name"), 
						"", 
						diff,
						diff,
						residue,
						position
						)
					);
			if( residue != null ){
				if( fixedAA[residueStr.charAt(0)-'A'] != 0 ) return 0;
				AminoAcid.modifiedAminoAcidMass( residueStr.charAt(0), diff );
				MSMass.modifiedAminoAcidMass( residueStr.charAt(0), diff );
				fixedAA[residueStr.charAt(0)-'A']= diff;
			}
		}
		return 1;
	}
	
	// since ver 1.0 Cashing Gap I
	public PTMSearchResult xxxsearchPTM( Sequence seq, double massDiff, PTMPosition position )
	{		
		PTMSearchResult searchResult;
			
		// Find result from hash table first
		searchResult = hashTable.get(new GapKey(seq, massDiff, position));
		if ( searchResult != null ){
			return searchResult;
		}

		// Hash miss
		ArrayList<PTMRun> newGapInterpret = new ArrayList<PTMRun>();
		
		// interpreted with no PTM
		// changed by timebird for mass correction

		if( Constants.fEqual(massDiff, 0.) )
		{
			searchResult = new PTMSearchResult(newGapInterpret, true);
			hashTable.put(new GapKey(seq, massDiff, position), searchResult);
			return searchResult;
		}

		// PTM search
		this.seq 		= seq;
		this.massDiff	= massDiff;
		this.position	= position;
		this.result		= newGapInterpret;
		findPTM_DFS( 0.0, 0, 0, 0 );
		this.result		= null;		// release reference from PTMDB class
		
		searchResult = new PTMSearchResult(newGapInterpret, newGapInterpret != null && newGapInterpret.size() > 0); 
		// add result to hash table
		hashTable.put(new GapKey(seq, massDiff, position), searchResult);
		
		return searchResult;
	}
	private double[][] modRange = new double[2][26];
	private double[][] ntermModRange = new double[2][26];
	private double[][] ctermModRange = new double[2][26];
	
	public double minimumModifiedMass( Sequence seq, int pos ){
		double min1=0, min2=0, extra = 0;
		for( AminoAcid aa : seq ){
			int aix= aa.getIndex();
			if( min1 > modRange[0][aix] ) {
				min2 = min1;
				min1 = modRange[0][aix];
			}
			else if( min2 > modRange[0][aix] ) min2 = modRange[0][aix];
		}
		
		if( pos == 1 ){
			int aix= seq.get(0).getIndex();
			if( min1 > ntermModRange[0][aix] ) {
				min2 = min1;
				min1 = ntermModRange[0][aix];
			}
			else if( min2 > ntermModRange[0][aix] ) min2 = ntermModRange[0][aix];
		}
		else if( pos == 2 ){
			int aix= seq.get(seq.size()-1).getIndex();
			if( min1 > ctermModRange[0][aix] ) {
				min2 = min1;
				min1 = ctermModRange[0][aix];
			}
			else if( min2 > ctermModRange[0][aix] ) min2 = ctermModRange[0][aix];
		}
		if( seq.size() == 1 ) return min1;
		return min1 + min2;
	}
	
	public double maximumModifiedMass(Sequence seq, int pos){
		double max1=0, max2=0;
		for( AminoAcid aa : seq ){
			int aix= aa.getIndex();
			if( max1 < modRange[1][aix] ) {
				max2 = max1;
				max1 = modRange[1][aix];
			}
			else if( max2 < modRange[1][aix] ) max2 = modRange[1][aix];
		}
		if( pos == 1 ){
			int aix= seq.get(0).getIndex();
			if( max1 < ntermModRange[1][aix] ) {
				max2 = max1;
				max1 = ntermModRange[1][aix];
			}
			else if( max2 < ntermModRange[1][aix] ) max2 = ntermModRange[1][aix];
		}
		else if( pos == 2 ){
			int aix= seq.get(seq.size()-1).getIndex();
			if( max1 < ctermModRange[1][aix] ) {
				max2 = max1;
				max1 = ctermModRange[1][aix];
			}
			else if( max2 < ctermModRange[1][aix] ) max2 = ctermModRange[1][aix];
		}
		if( seq.size() == 1 ) return max1;
		return max1 + max2;
	}
	
	public void setPTMDiagnosticIon() {		
		for (PTM ptm : this) {
			
			//setting common modifications
			if( Math.abs(15.994915-ptm.getMassDifference()) < 0.01 ){ 
				if( ptm.getAbbAA() =='M' ) {
					ptm.setPenalty(0.5);
					ptm.setModCount(0);
					AminoAcid.canBeEasilyModified('M');
					ptm.setNeutralLoss( ptm.getMassDifference()+48.003371 );
				}
				if( ptm.getAbbAA() =='W' ) ptm.setPenalty(0.8);	
			}
			
			if( ptm.getAbbAA() =='Q' && ptm.getPTMPosition() == PTMPosition.ANY_N_TERM && Constants.round(ptm.getMassDifference()+Constants.NTERM_FIX_MOD) == -17 ){
				ptm.setPenalty(0.5);
				ptm.setModCount(0);
			}

			//setting isobaric labeling 
			if( Constants.reporterMassOfIsobaricTag != null && Math.abs(Constants.reporterMassOfIsobaricTag[0]-ptm.getMassDifference()) < 0.01 ){
				if( ptm.getAbbAA() == 'S' || ptm.getAbbAA() == 'T' ) {
					if( "itraq4plex".compareTo(Constants.isobaricTag) == 0 ) ptm.setDiagnosticIon(163.1199);
					ptm.setPenalty(0.5);
					ptm.setModCount(0);
					AminoAcid.canBeEasilyModified(ptm.getAbbAA());
				}
			}
			
			//setting specific modifications
			if( ptm.getName().startsWith("Acetyl") ) {
				if( ptm.getAbbAA() =='K' ) ptm.setDiagnosticIon(126.0913);
				if( ptm.getPTMPosition() == PTMPosition.PROTEIN_N_TERM ) ptm.setPenalty(0.5);				
				
				if( "Acetyl".compareToIgnoreCase(Constants.enrichedModification) == 0 ) {
					ptm.setPenalty(0.5);					
					if( ptm.getAbbAA() =='K' ) {
						ptm.setModCount(0);
						AminoAcid.canBeEasilyModified('K');
					}
					if( ptm.getPTMPosition() == PTMPosition.PROTEIN_N_TERM ) {
						ptm.setModCount(0);
					}
				}
			}
			
			if( "Phospho".compareTo(ptm.getName()) == 0 ){
				if( ptm.getAbbAA() =='S' || ptm.getAbbAA() =='T' ) {
					ptm.setNeutralLoss(ptm.getMassDifference()+Constants.H2O);
					if( "Phospho".compareToIgnoreCase(Constants.enrichedModification) == 0 ) {
						ptm.setPenalty(0.5);
						ptm.setModCount(0);
						AminoAcid.canBeEasilyModified(ptm.getAbbAA());
					}
				}
			}
		}
	}
	
}


