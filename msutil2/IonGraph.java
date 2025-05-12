package msutil;

import java.util.ArrayList;

import moda.ThreadLocalMutables;
import modi.Constants;
import modi.Mutables;
import modi.PTM;

@SuppressWarnings("unused")
public abstract class IonGraph extends ArrayList<IonNode> {
	
	static final int NTermType = 0;
    static final int CTermType = 1;
	static final int SinglyCharged = 1;
    static final int DoublyCharged = 2;
    static final int TriplyCharged = 3;
	static final double supportingFactor = 0.9;
	
	final String 	sequence;
	final double[] ptmMass;
	final PTM[] 	 ptmList;
	final int charge;
	
	double 	calculatedMW = 0;	
	double 	modifiedResd = 0;
	
	int 	rankScore = 0;
	int 	coverScore =0;
	double	prob = 0;
	double	massError = 0;
	
	double ionCoverage = 0;
	double seqCoverage = 0;
	
	int	   maxConsecutive = 0; 	
	int	   maxAAs = 0; //AA Length, not inos

	public double getCalculatedMW(){ return calculatedMW; }
	public double getMassError(){ return massError; }
	public int getRankScore(){ return rankScore; }
	public int getCoverScore(){ return coverScore; }
	public int getModifiedResd(){ return Constants.round(modifiedResd); }
	public double getProb(){ return prob; }
	public double getIonCoverage(){ return ionCoverage; }
	public double getSeqCoverage(){ return seqCoverage; }
	public int getMaxConsecutiveAALength(){ return maxAAs; }
	
	public IonGraph( String peptide, double[] ptms, PGraph graph ){	
		sequence = peptide;
		ptmMass  = ptms;
		charge   = graph.getCharge();
		ptmList = null;
		construct();
	}
	
	public IonGraph( String peptide, double[] ptms, PTM[] pList, PGraph graph ){	
		sequence = peptide;
		ptmMass = ptms;
		charge = graph.getCharge();
		ptmList = pList;
		construct();
	}
	
	protected abstract void construct();
	public abstract void setScore( PGraph graph );
	protected abstract void setSecondaryScore( PGraph graph );
	public abstract void evaluateMatchQuality( PGraph graph );
	
	public void setMatchCoverage(){

		int len = sequence.length()-1;
		int[][] ionHits = new int[2][len];
		int[] seqHits = new int[len];
		
		for( IonNode node : this ){
			ionHits[node.type][node.index]++;
			seqHits[node.index]++;
		}
		
		int bbHit=0, yyHit=0, seqHit=0;
		for( int i=0 ; i<len; i++ ){
			if( ionHits[NTermType][i] != 0 ) bbHit++; 
			if( ionHits[CTermType][i] != 0 ) yyHit++;
			if( seqHits[i] != 0 ) seqHit++;
		}
		if( bbHit < yyHit ){
			bbHit = yyHit;
		}
		
		int maxConsec = 0, localsec = 0;
		for( int i=0 ; i<len; i++ ){// BION
			if( ionHits[NTermType][i] != 0 ) localsec++;
			else{
				if( maxConsec < localsec ) maxConsec = localsec;
				localsec = 0;		
			}
		}
		if( maxConsec < localsec ) maxConsec = localsec;
		
		localsec = 0;
		for( int i=0 ; i<len; i++ ){// YION
			if( ionHits[CTermType][i] != 0 ) localsec++;
			else{
				if( maxConsec < localsec ) maxConsec = localsec;
				localsec = 0;
			}
		}
		if( maxConsec < localsec ) maxConsec = localsec;
		
		maxConsecutive = maxConsec;
		maxAAs  = getMaxAAs(ionHits);
		
		double modifiedMass = 0;
		for(int i=0; i<ptmMass.length ; i++){
			if( ptmMass[i] != 0 ){
				modifiedMass += ptmMass[i];
			}
		}
		
		int modAcid = ( modifiedMass < (ThreadLocalMutables.get().precursorTolerance) )?  0 : (int)Math.ceil( modifiedMass / 110 );
		double penalty = modAcid + modifiedResd/2;
		seqCoverage = (double)seqHit/(len+penalty);		
		ionCoverage = (double)bbHit /(len+penalty);
	}
	
	protected int getMaxAAs(int[][] ionHits){

		int whichIonIntense = 2; //0:bion, 1:yion, 2:dont know	

		int len = ionHits[0].length;
		int maxConsec = 0, localsec = 0;		
		for( int i=0 ; i<len; i++ ){// YION
			if( ionHits[CTermType][i] != 0 ) localsec++;
			else{
				if( maxConsec < localsec ) maxConsec = localsec;
				localsec = -1;
			}
		}
		localsec++;
		if( maxConsec < localsec ) maxConsec = localsec;
		if( whichIonIntense == CTermType ) return maxConsec;

		if( whichIonIntense == NTermType ) maxConsec = 0;
		localsec = 0;	
		for( int i=0 ; i<len; i++ ){// BION
			if( ionHits[NTermType][i] != 0 ) localsec++;
			else{
				if( maxConsec < localsec ) maxConsec = localsec;
				localsec = -1;		
			}
		}
		localsec++;
		if( maxConsec < localsec ) maxConsec = localsec;		

		return maxConsec;
	}
	
	protected void reCalibrateMZ(){
		if( this.size() < 2 ) return; 
		
		double[] mzError = regression();	
		massError = calcStdDeviation( mzError );		
		double cutOutlier = massError * 4;
		
		int mis = 0;
		for(int i=0; i<mzError.length; i++){
			if( Math.abs(mzError[i]) > cutOutlier ){
				this.remove(i-mis);
				mis++;
			}		
		}
		if( mis > 0 ){
			double[] mzError2 = regression();		
			massError = calcStdDeviation( mzError2 );	
		}	  
	}
	
	protected double[] regression(){
		int i;
		
		double[] mzError = new double[this.size()];
		double aConst, bSlope=1;// regression ef.
		double x_cap=0, y_cap=0, temp=0;
		
		i=0; ////////////////////////////////
		while( i < this.size() ){		
			x_cap += ( this.get(i).intensity * this.get(i).mz );
			y_cap += ( this.get(i).intensity * this.get(i).observed );
			temp += this.get(i).intensity;		
			i++;		
		}
		if( temp == 0 ) return mzError;		
		x_cap /= temp;
		y_cap /= temp;

		i=0; ////////////////////////////////
		double divSum=0;
		while( i < this.size() ){
			bSlope += ( this.get(i).intensity * (this.get(i).mz-x_cap) * (this.get(i).observed-y_cap) );
			divSum += ( this.get(i).intensity * Math.pow((this.get(i).mz-x_cap), 2) );
			i++;
		}
		
		bSlope= ( bSlope == 0 )? 1 : bSlope/divSum;
		aConst= y_cap - ( bSlope*x_cap );
		
		i=0; ////////////////////////////////
		while( i < this.size() ){
			mzError[i] = this.get(i).observed - ( aConst + bSlope*this.get(i).mz );
			i++;
		}

		return mzError;	
	}
	
	protected double calcStdDeviation(double[] mzError){
		double sd= 0;
		for(int i=0; i<mzError.length; i++)
			sd +=  Math.pow(mzError[i], 2);
		sd = sd / (mzError.length-1);
		return Math.sqrt(sd);	
	}
	
}













