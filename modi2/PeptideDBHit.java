package modi;

import java.io.Serializable;

public class PeptideDBHit implements Serializable {
	private final Peptide 		hitPeptide;
	private int 			start;
	private int 			end;
	private final double			nTermOffset;
	private final double			cTermOffset;
	
	public int			getStart()			{ return start; }
	public int			getEnd()			{ return end; }
	public double		getNTermOffset()	{ return nTermOffset; }
	public double		getCTermOffset()	{ return cTermOffset; }
	public Peptide		getHitPeptide()		{ return hitPeptide; }
	
	public void			setStart(int start)	{ this.start = start; }
	public void			setEnd(int end)		{ this.end = end; }
	@SuppressWarnings("BooleanMethodIsAlwaysInverted")
    public boolean		isNtermIon()		{ return start == 0; }
	@SuppressWarnings("BooleanMethodIsAlwaysInverted")
    public boolean		isCtermIon()		{ return end == hitPeptide.size()-1; }
	
	public PeptideDBHit( Peptide hitPeptide, int start, int end, double nTermOffset, double cTermOffset )
	{
		this.hitPeptide		= hitPeptide;
		this.start			= start;
		this.end			= end;
		this.nTermOffset	= nTermOffset;
		this.cTermOffset	= cTermOffset;		
	}

}
