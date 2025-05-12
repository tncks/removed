package modi;

class DiagnosticPTMPeak extends Peak
{
	final double		prob;
	public DiagnosticPTMPeak(double mass, double prob)
	{
		super(-1, mass, 0, 1, PeakProperty.DIAGNOSTIC_PTM_ION);
		this.prob = prob;
	}
}

