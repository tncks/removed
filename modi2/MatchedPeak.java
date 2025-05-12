package modi;

public class MatchedPeak extends Peak {
	final AminoAcid NTermResidue;
	final AminoAcid CTermResidue;
	
	public MatchedPeak(int index, double mass, double intensity, int charge, PeakProperty property, 
			AminoAcid NTermResidue, AminoAcid CTermResidue)
	{
		super(index, mass, intensity, charge, property);
		this.NTermResidue = NTermResidue;
		this.CTermResidue = CTermResidue;
	}
	
	public AminoAcid getNTermResidue()	{ return NTermResidue; }
	public AminoAcid getCTermResidue()	{ return CTermResidue; }
	
}
