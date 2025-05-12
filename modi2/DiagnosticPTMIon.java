package modi;


public class DiagnosticPTMIon {
	final String		name;
	final AminoAcid	residue;
	final double		mass;
	final double		prob;
	
	public DiagnosticPTMIon(String name, AminoAcid residue, double mass, double prob)
	{
		this.name = name;
		this.residue = residue;
		this.mass = mass;
		this.prob = prob;
	}

	public	double		getMass()		{ return mass; }
	public	double		getProb()		{ return prob; }
}

