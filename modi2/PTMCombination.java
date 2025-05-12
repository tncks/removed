package modi;

public class PTMCombination {
	String 		ptmComb = "";
	final double[] 	ptms;
	final PTM[] 		ptmList;
	
	public PTMCombination(int len){
		ptms = new double[len];
		ptmList = new PTM[len];
	}
	
	public int hashCode(){ 
		return ptmComb.hashCode();
	}
	
	public boolean equals(Object o){		
		if(!(o instanceof PTMCombination p))
			   return false;
        return this.ptmComb.equals(p.ptmComb);
    }
}
