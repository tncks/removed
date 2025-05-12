package modi;



public class ScanCap implements Comparable<ScanCap> {
	private final String 	title;
	private final double 	pmz;
	private final double 	neutralMW;
	private final int 	charge;
	private int 	scanNo;
	private long 	offset;
	
	public ScanCap(String title, double pmz, int charge){
	
		this.title 		= title;	
		this.pmz 		= pmz;
		this.charge 	= charge;
		this.neutralMW 	= (pmz - Constants.Proton)*charge;
	}
	
	public ScanCap(String title, int sn, double pmz, int charge){		
		this.title 		= title;	
		this.scanNo		= sn;
		this.pmz 		= pmz;
		this.charge 	= charge;
		this.neutralMW 	= (pmz - Constants.Proton)*charge;
	}

	public void setOffset(long offset){ this.offset = offset; }
	public String getTitle(){ return title; }
	public int getScanNumber(){ return scanNo; }
	public double getObservedMW(){ return neutralMW; }
	public double getPMZ(){ return pmz; }
	public int getCharge(){ return charge; }
	public long getOffset(){ return offset; }
	

	
	private static class RawP implements Comparable<RawP> {
		final double mz;
		final double it;
		public RawP(double m, double i){
			mz=m;
			it=i;
		}	
		public int compareTo(RawP p) {
			if( mz > p.mz ) return 1;
			else if( mz < p.mz ) return -1;
			else return 0;
		}
	}
	
	public int compareTo(ScanCap s) 
	{
		if( this.neutralMW > s.neutralMW ) return 1;
		else if( this.neutralMW < s.neutralMW ) return -1;
		
		if( this.charge > s.charge ) return 1;
		else if( this.charge < s.charge ) return -1;
		else return 0;
	}
	

}
