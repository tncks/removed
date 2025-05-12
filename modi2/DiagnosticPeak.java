package modi;

class DiagnosticPeak extends Peak
{
	final double		probAAOverIon;
	final double		probIonOverAA;
	public DiagnosticPeak(double mass, double probAAOverIon, double probIonOverAA)
	{
		super(-1, mass, 0, 1, PeakProperty.DIAGNOSTIC_ION);
		this.probAAOverIon = probAAOverIon;
		this.probIonOverAA = probIonOverAA;
	}
}

