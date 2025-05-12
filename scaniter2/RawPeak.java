package scaniter;

public class RawPeak implements Comparable<RawPeak> {

	final double mz;
	final double it;
	
	public RawPeak(double m, double i){
		mz=m;
		it=i;
	}
	
	public int compareTo(RawPeak p) {
		if( mz > p.mz ) return 1;
		else if( mz < p.mz ) return -1;
		else return 0;
	}
}


