package msutil;

import java.util.Comparator;

public class TNode implements Comparable<TNode> {
	final char 	ionType; // B, Y, A, C, X, Z
	final int 	charge;
	final int 	rank;
	final double 	theoMZ;
    final double obsvMZ;
	final double 	obsvIT;
	double 	error;
	
	public TNode(char type, int cs, double tmz, double omz, double oit, int r){
		ionType = type;
		charge = cs;
		theoMZ = tmz;
		obsvMZ = omz;
		obsvIT = oit;
		rank = r;
	}
	
	public int compareTo(TNode o) {
		if( this.theoMZ > o.theoMZ )
			return 1;
		else if( this.theoMZ == o.theoMZ )
			return 0;
		else
			return -1;
	}	
	
}
class TNodeRankComparator implements Comparator<TNode>{
	public int compare(TNode x1, TNode x2){
		if( x1.rank > x2.rank ) return 1;
		else if( x1.rank < x2.rank ) return -1;
		else return 0;
	}	
	public boolean equals(TNode x1, TNode x2){
		return x1 == x2;
	}
}