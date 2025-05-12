package msutil;

import modi.Constants;

@SuppressWarnings("unused")
public class IonNode implements Comparable<IonNode> {
	final int 	index;
	final int 	type; // nterm=1, cterm=2
	final int 	charge;
	final double 	mz; // theoretical
	final int 	priority;
	double 	observed;
	double 	intensity;
	
	public IonNode(int i, int t, int cs, double m, int p){
		index = i;
		type = t;
		mz = m;
		charge =cs;
		priority = p;
	}
	
	public int getIndex(){ return index; }
	public int getCharge(){ return charge; }
	public int getType(){ return type; }
	public double getMZ(){ return mz; }
	public int getPriority(){ return priority; }
	public double getObserved(){ return observed; }
	public double getMZError(){ return Math.abs( observed - mz ); }
	public double getMW(){ 
		return (mz-Constants.Proton)*charge; 
	}
	
	
	public int compareTo(IonNode o) {
		if( this.priority > o.priority ) return 1;
		else if( this.priority == o.priority ){
			if( this.mz > o.mz ) return 1;
			else if( this.mz < o.mz ) return -1;
			else return 0;
		}
		else return -1;
	}	
}
