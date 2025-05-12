package modi;

import java.util.ArrayList;

public class AminoAcid implements Comparable<AminoAcid> {

	private char 	residue;
	private double	monoMass;
	private double	avgMass;
	private String	fullName;
	private int		index;		// used for Bucket indexing
	private boolean isLabelled = false;
	
	private static int indexSize = 0;
	
	public int compareTo( AminoAcid o ){
		return this.residue - o.residue;
	}
	
	public AminoAcid(char residue, double monoMass, double avgMass, String fullName, int index )
	{
		this.residue	= residue;
		this.monoMass	= monoMass;
		this.avgMass	= avgMass;
		this.fullName	= fullName;
		this.index		= index;
		
		if (AminoAcid.indexSize<index+1) AminoAcid.indexSize = index+1;
	}
	
	public String toString()				
	{
		if(index == 24)
			return "N'";
		else if(index == 25)
			return "Q'";
		return String.valueOf(residue); 
	}

	public char			getResidue()		{ return residue; }
	public double		getMass()			{ return monoMass; }
	public double		getMonoMass()		{ return monoMass; }
	public double		getAvgMass()		{ return avgMass; }
	public String		getFullName()		{ return fullName; }
	public int			getIndex()			{ return index; }
	public boolean		isLabelled()		{ return isLabelled; }
	public static int	getIndexSize()		{ return indexSize; }
	public static double getMass(char residue) { return getAminoAcid(residue).getMass(); }
	
	public static ArrayList<AminoAcid> getCode(double mass, double massTolerance)
	{
		ArrayList<AminoAcid> resultCode = new ArrayList<AminoAcid>();
		// linear search, might be optimised later 
		for(int i=0; i<aaTable.length; i++)  {
			if( aaTable[i] == null )
				continue;
	//		if( aaTable[i].residue == 'L' || aaTable[i].residue == 'Q' )
	//			continue;
			
			else if( aaTable[i].getMass() < mass + massTolerance &&
					aaTable[i].getMass() > mass - massTolerance )
				resultCode.add( aaTable[i] );
		}
		return resultCode;
	}

	public static double getMaxMonoMass() { return maxAAMass; }
	public static double getMinMonoMass() { return minAAMass; }
	
	public static AminoAcid getAminoAcid( char residue ) {
		if (residue>='a' && residue<='z')
			return aaTable[residue-'a'];
		else if (residue>='A' && residue<='Z')
			return aaTable[residue-'A'];
		
		return null;
	}
	
	public	static AminoAcid getAminoAcid(int index) {
		for(AminoAcid aa : aaTable)
			if(aa != null && aa.getIndex() == index)
				return aa;

		return null;
	}
	
	private static double minAAMass = 57.02146;
	private static double maxAAMass = 186.07931;
	// Static table containing Predefined Amino Acids	
	private static AminoAcid [] aaTable = 
	{		
		new AminoAcid('A',	71.037114,	71.08,	"alanine",			0),
		null, //new AminoAcid('B',	0,			0,		"_invalid_" ),
	//	new AminoAcid('C',	103.00919+57.0215,	103.1+57.0213,	"cysteine",			1),
		new AminoAcid('C',	103.00919,	103.1,	"cysteine",			1),
		new AminoAcid('D',	115.02694,	115.1,	"aspartate",		2),
		new AminoAcid('E',	129.04259,	129.1,	"glutamate",		3),
		new AminoAcid('F',	147.06841,	147.2,	"phenylalanine", 	4),
		new AminoAcid('G',	57.021464,	57.05,	"glycine", 			5),
		new AminoAcid('H',	137.05891,	137.1,	"histidine", 		6),
		new AminoAcid('I',	113.08406,	113.2,	"isoleucine", 		7),
		null,
		new AminoAcid('K',	128.09496,	128.2,	"lysine", 			8),
		new AminoAcid('L',	113.08406,	113.2,	"leucine", 			9),
		new AminoAcid('M',	131.04048,	131.2,	"methionine", 		10),
		new AminoAcid('N',	114.04293,	114.1,	"asparagine", 		11),
		null, //new AminoAcid('O',	0,			0,		"_invalid_"),
		new AminoAcid('P',	97.052764,	97.12,	"proline", 			12),
		new AminoAcid('Q',	128.05858,	128.1,	"glutamine", 		13),
		new AminoAcid('R',	156.10111,	156.2,	"arginine", 		14),
		new AminoAcid('S',	87.032029,	87.08,	"serine", 			15),
		new AminoAcid('T',	101.04768,	101.1,	"threonine", 		16),
		null, //new AminoAcid('U',	151.00919,	151.1,	"selenocysteine",	20),
		new AminoAcid('V',	99.068414,	99.07,	"valine",			17),
		new AminoAcid('W',	186.07931,	186.2,	"tryptophan",		18),
		null, //new AminoAcid('X',	0,			0,		"_invalid_"),
		new AminoAcid('Y',	163.06333,	163.2,	"tyrosine",			19),
		null, //new AminoAcid('Z',	0,			0,		"_invalid_"),
		null, //new AminoAcid('*',	0,			0,		"_invalid_"),
		null //new AminoAcid('-',	0,			0,		"_invalid_"),
	//	new AminoAcid('@',	113.08406,	113.2, "leusine or isoleucine",	20),
	//	new AminoAcid('*',			0,		0, "N-terminal",			21),
	//	new AminoAcid('-',			0,		0, "C-terminal",			22),
	//	new AminoAcid('&',			0,		0, "Modified Amino Acid", 	23),
	//	new AminoAcid('~',	115.04293,	115.1,	"deamidated asparagine",24),
	//	new AminoAcid('`',	129.05858,	129.1,	"deamidated glutamine", 25),
	};
	
	public static void modifiedAminoAcidMass(char AA, double fixedModification){
		aaTable[AA-'A'].monoMass += fixedModification;	
		aaTable[AA-'A'].avgMass += fixedModification;
		if( aaTable[AA-'A'].monoMass > maxAAMass ) maxAAMass = aaTable[AA-'A'].monoMass;
		if( aaTable[AA-'A'].monoMass < minAAMass ) minAAMass = aaTable[AA-'A'].monoMass;
		if( AA != 'C' ) aaTable[AA-'A'].isLabelled = true;
	}
	public static void canBeEasilyModified(char AA){
		aaTable[AA-'A'].isLabelled = true;
	}
	
	public static char getAminoAcidBy3Letter( String threeL ) {
		if( threeL.compareTo("Ala") == 0 ) return 'A';
		if( threeL.compareTo("Cys") == 0 ) return 'C';
		if( threeL.compareTo("Asp") == 0 ) return 'D';
		if( threeL.compareTo("Glu") == 0 ) return 'E';
		if( threeL.compareTo("Phe") == 0 ) return 'F';
		if( threeL.compareTo("Gly") == 0 ) return 'G';
		if( threeL.compareTo("His") == 0 ) return 'H';
		if( threeL.compareTo("Xle") == 0 ) return 'I';
		if( threeL.compareTo("Lys") == 0 ) return 'K';
	//	if( threeL.compareTo("Xle") == 0 ) return 'L';
		if( threeL.compareTo("Met") == 0 ) return 'M';
		if( threeL.compareTo("Asn") == 0 ) return 'N';
		if( threeL.compareTo("Pro") == 0 ) return 'P';
		if( threeL.compareTo("Gln") == 0 ) return 'Q';
		if( threeL.compareTo("Arg") == 0 ) return 'R';
		if( threeL.compareTo("Ser") == 0 ) return 'S';
		if( threeL.compareTo("Thr") == 0 ) return 'T';
		if( threeL.compareTo("Val") == 0 ) return 'V';
		if( threeL.compareTo("Trp") == 0 ) return 'W';
		if( threeL.compareTo("Tyr") == 0 ) return 'Y';
		return 'X';
	}
	
	public static String toRRIX(){
		StringBuffer des = new StringBuffer(indexSize+"|");
		for(int i=0; i<aaTable.length; i++) {
			if( aaTable[i] == null )
				continue;	
			des.append( aaTable[i].toString()+String.format("|%.4f|", aaTable[i].getMass()) );
		}		
		return des.toString();
	}
}
