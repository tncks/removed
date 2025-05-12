package modi;

@SuppressWarnings("ClassCanBeRecord")
class DiagnosticIon {
	final AminoAcid	aa;
	final double		mass;
	final double		probAAOverIon;
	final double		probIonOverAA;
	
	public DiagnosticIon(AminoAcid aa, double mass, double probAAOverIon, double probIonOverAA)
	{
		this.aa = aa;
		this.mass = mass;
		this.probAAOverIon = probAAOverIon;
		this.probIonOverAA = probIonOverAA;
	}
	public	AminoAcid	getAminoAcid()		{ return aa; }
	public	double		getMass()			{ return mass; }
	public	double		getProbAAOverIon()	{ return probAAOverIon; }
	public	double		getProbIonOverAA()	{ return probIonOverAA; }
}