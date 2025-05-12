package msutil;

import java.util.ArrayList;
import java.util.Collections;

import modi.Constants;

public class TGraph extends ArrayList<TNode> {
	
	double 	std;
	int 	multiMatches = 0;
	double		modified = 0;
	final double 	calculatedMW;
	
	public double getCalculatedMW(){ return calculatedMW; }
	
	public TGraph( String peptide, double[] ptms, PGraph graph ){
		int len = peptide.length()-1;
		int i;
		
		double ions= Constants.B_ION_OFFSET;
		for(i=0; i<len; i++ ){
			ions += MSMass.getAAMass( peptide.charAt(i) );
			ions += ptms[i];
			if( ptms[i] != 0 ) modified += ptms[i];
			PNode p = graph.getMatchedPNode(ions);
			if( p != null ){
				this.add( new TNode('B', 1, ions, p.mass, p.norm, p.localRank) );
			}
			if( i == 1 ){
				PNode ap = graph.getMatchedPNode(ions- Constants.A_ION_OFFSET);
				if( ap != null ){
					this.add( new TNode('A', 1, ions- Constants.A_ION_OFFSET, ap.mass, ap.norm, ap.localRank) );
				}
			}
		}
		ions += MSMass.getAAMass( peptide.charAt(i) );
		ions += ptms[i];
		if( ptms[i] != 0 ) modified += ptms[i];
		calculatedMW = ions - Constants.B_ION_OFFSET + Constants.H2O;
		
		ions= Constants.Y_ION_OFFSET;
		for(i=len; i>0; i-- ){
			ions += MSMass.getAAMass( peptide.charAt(i) );
			ions += ptms[i];
			PNode p = graph.getMatchedPNode(ions);
			if( p != null ){
				this.add( new TNode('Y', 1, ions, p.mass, p.norm, p.localRank) );
			}
		}
		Collections.sort(this);
	}
	public TGraph( String peptide, int[] ptms, PGraph graph ){
		int len = peptide.length()-1;
		int i;
		
		double ions= Constants.B_ION_OFFSET;
		for(i=0; i<len; i++ ){
			ions += MSMass.getAAMass( peptide.charAt(i) );
			ions += ptms[i];
			PNode p = graph.getMatchedPNode(ions);
			if( p != null ){
				this.add( new TNode('B', 1, ions, p.mass, p.norm, p.localRank) );
			}
			if( i == 1 ){
				PNode ap = graph.getMatchedPNode(ions- Constants.A_ION_OFFSET);
				if( ap != null ){
					this.add( new TNode('A', 1, ions- Constants.A_ION_OFFSET, ap.mass, ap.norm, ap.localRank) );
				}
			}
		}
		
		ions += MSMass.getAAMass( peptide.charAt(i) );
		ions += ptms[i];
		
		calculatedMW = ions - Constants.B_ION_OFFSET + Constants.H2O;
		
		ions= Constants.Y_ION_OFFSET;
		for(i=len; i>0; i-- ){
			ions += MSMass.getAAMass( peptide.charAt(i) );
			ions += ptms[i];
			PNode p = graph.getMatchedPNode(ions);
			if( p != null ){
				this.add( new TNode('Y', 1, ions, p.mass, p.norm, p.localRank) );
			}
		}
		Collections.sort(this);
	}
	public TGraph( String peptide, PGraph graph ){
		int len = peptide.length()-1;
		int i;
		
		double ions= Constants.B_ION_OFFSET;
		for(i=0; i<len; i++ ){
			ions += MSMass.getAAMass( peptide.charAt(i) );
			PNode p = graph.getMatchedPNode(ions);
			if( p != null ){
				this.add( new TNode('B', 1, ions, p.mass, p.norm, p.localRank) );
			}
			if( i == 1 ){
				PNode ap = graph.getMatchedPNode(ions- Constants.A_ION_OFFSET);
				if( ap != null ){
					this.add( new TNode('A', 1, ions- Constants.A_ION_OFFSET, ap.mass, ap.norm, ap.localRank) );
				}
			}
		}
		
		ions += MSMass.getAAMass( peptide.charAt(i) );
		
		calculatedMW = ions - Constants.B_ION_OFFSET + Constants.H2O;
		
		ions= Constants.Y_ION_OFFSET;
		for(i=len; i>0; i-- ){
			ions += MSMass.getAAMass( peptide.charAt(i) );
			PNode p = graph.getMatchedPNode(ions);
			if( p != null ){
				this.add( new TNode('Y', 1, ions, p.mass, p.norm, p.localRank) );
			}
		}
		Collections.sort(this);
	}
	
	public double getExplainedScore( PGraph graph ){ // iontype(n/c), charge(1/2)
		graph.refresh();
		ArrayList<Integer> mats = new ArrayList<>();
		
		double score = 0;						
		for( TNode tn : this ){
			int index = graph.getIndexOfMatchedPeak( tn.theoMZ );
			if( index != -1 ){
				score += tn.obsvIT;
				mats.add( index );
			}
		}
		
		for(int i=0; i<mats.size(); i++){
			score += graph.getSupportingScore( mats.get(i), 1 );
		}		
		return score;
	}
	
	
	public double getRankScore( PGraph graph ){
		graph.refresh();
		ArrayList<Integer> mats = new ArrayList<>();
		double score = 0;						
		for( TNode tn : this ){
			int index = graph.getIndexOfMatchedPeak( tn.theoMZ );
			if( index != -1 ){
				score += tn.obsvIT;
				mats.add(index);
			}
		}
		for(int i=0; i<mats.size(); i++){
			score += graph.getSupportingScore( mats.get(i), 1 );
		}		
		return score;
	}
	
	public void process(){
		reomveOutlier();
		reomveMultiMatches();
	}
	
	private void reomveMultiMatches(){	
		int rev = 0;
		for(int i=0; i<this.size()-1; i++){					
			if( this.get(i).obsvMZ == this.get(i+1).obsvMZ ){
				if( this.get(i).error >  this.get(i+1).error ){
					this.remove(i);
					i--;
				}
				else this.remove(i+1);
				rev++;	
			}
		}	
		multiMatches = rev;
	}
	
	private void reomveOutlier(){
		calibrateMZ();
		calcStdDeviation();
		
		double cutOutlier = this.std * 3;
		for(int i=0; i<this.size(); i++){
			if( Math.abs(this.get(i).error) >  cutOutlier ){
				this.remove(i);
				i--;
			}			
		}
		calcStdDeviation();
	}
	
	private int calibrateMZ(){
		int i;
		double aConst, bSlope=0;// regression ���.
		double x_cap=0, y_cap=0, temp=0;
		
		i=0; ////////////////////////////////
		while( i < this.size() ){		
			x_cap += ( this.get(i).obsvIT * this.get(i).theoMZ );
			y_cap += ( this.get(i).obsvIT * this.get(i).obsvMZ );
			temp += this.get(i).obsvIT;
			i++;		
		}
		if( temp == 0 ) return 0;		
		x_cap /= temp;
		y_cap /= temp;

		i=0; ////////////////////////////////
		double divSum=0;
		while( i < this.size() ){
			bSlope += ( this.get(i).obsvIT * (this.get(i).theoMZ-x_cap) * (this.get(i).obsvMZ-y_cap) );
			divSum += ( this.get(i).obsvIT * Math.pow((this.get(i).theoMZ-x_cap), 2) );
			i++;
		}
		
		bSlope= ( bSlope == 0 )? 1 : bSlope/divSum;
		aConst= y_cap - ( bSlope*x_cap );
		
		int maxR=0;
		i=0; ////////////////////////////////
		temp=0;
		while( i < this.size() ){
			this.get(i).error = this.get(i).obsvMZ - ( aConst + bSlope*this.get(i).theoMZ ) ;
			if( this.get(i).error > temp ){
				maxR= i;
				temp = this.get(i).error;
			}
			i++;
		}
		return maxR;	
	}
	
	private void calcStdDeviation(){
		if( this.size() < 2 ) { 
			this.std=0; 
			return; 
		}
		int i=0;
		double sd= 0;

		while( i < this.size() ){		
			sd +=  Math.pow(this.get(i).error, 2);
			i++;
		}
		sd = sd / (this.size()-1);
		this.std = Math.sqrt(sd);	
	}
	
	public void sortByRank(){
		this.sort(new TNodeRankComparator());
	}
	
}
