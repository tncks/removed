package msutil;
import java.util.ArrayList;
import java.util.Collections;

import moda.ThreadLocalMutables;
import modi.Constants;
import modi.Mutables;

@SuppressWarnings("unused")
public class PGraph extends ArrayList<PNode>{
	private final double 	PMZ;
    private final double obsvMW;
    private double corrMW;
	private double	TIC = 0.;
	private final int 	charge;
	private static final int isotopeRatio = 5;

	public PGraph( double m, int c ){
		corrMW = obsvMW = m; 
		charge = c;
		PMZ = (obsvMW + charge*Constants.Proton)/charge;
	}
	
	public int getCharge(){ return charge; }
	public double getObservedMW(){ return obsvMW; }	
	public double getCorrectedMW(){ return corrMW; }
	
	public void ready( double x ){ 		
		TIC = x; 
		rankingArrange();
		Collections.sort(this);
	}
	public void refresh(){
		for( PNode x : this )
			x.assign( false );
	}
	
	public void setSelectedTIC(){ 
		TIC = 0;
		for( int i=0; i<this.size(); i++){			
			TIC += this.get(i).intensity;
		}
	}
	
	public double correctMW( boolean dynamicCorrection ){ 
		
		if( Constants.INSTRUMENT_TYPE == Constants.msms_type.QTOF ) return obsvMW;				
		if( (ThreadLocalMutables.get().precursorTolerance) <= Mutables.fragmentTolerance ) return obsvMW;
		
		PRM prmTable;
		if( charge > 2 ) prmTable= new PRMforHighCharge(this);
		else prmTable= new PRM(this);//*/
		
		double error = prmTable.massCorrection( dynamicCorrection );
		corrMW = obsvMW - error;			
		return corrMW; 
	}
	
	public void rankingArrange(){
		double preCount = -100.;
		int preRank = 100;
		for( int i=this.size()-1; i>-1; i-- ){
			
			int originRank = this.get(i).localRank;
			
			if( this.get(i).intensity == preCount ){
				this.get(i).localRank = preRank;
				this.get(i).norm = Constants.rNorm[preRank];
			}
			else{
				preCount = this.get(i).intensity;
				preRank = this.get(i).localRank;
			}
			
			if( originRank == 1 ){
				preCount = -100.;
				preRank = 100;
			}
		}
	}

	public double getExplainedTIC(){ 
		double mat=0;
		for( int i=0; i<this.size(); i++){
			if( this.get(i).assigned ) {
				mat += this.get(i).intensity;
			}
		}
		return mat/TIC;
	 }
	
	public double getPenaltySCoreForUnassignedIons(int N){ 		
		this.sortRevByPNodeIntensity();			
		double unmat=0;
		int pick = 0;
		for( int i=0; i < this.size() ; i++){		
			if( this.get(i).mass < PMZ && PMZ - this.get(i).mass < 50 ) continue;
			if( !this.get(i).assigned ) {
				unmat += this.get(i).norm;
			}		
			if( ++pick == N ) break;
		}
		unmat += N - pick;
		Collections.sort( this );
		return unmat;
	}
	public int getUnassignedTopN(int top){ 

		this.sortRevByPNodeIntensity();			
		int unmat=0, extra=0;
		for( int i=0; i<this.size() && i<top+extra ; i++){
			if( this.get(i).mass < PMZ && PMZ - this.get(i).mass < 50 ) {
				extra++;
				continue;			
			}
			if( !this.get(i).assigned ) {
				unmat++;
			}
		}//*/

		Collections.sort( this );
		return unmat;
	 }
	public double getUnassignedNorm(int lr){ 					
		
		double unmat=0;
		for( int i=0; i < this.size() ; i++){			
			if( this.get(i).mass < PMZ && PMZ - this.get(i).mass < 50 ) continue;			
			if( !this.get(i).assigned && this.get(i).localRank <= lr ) {
				unmat += this.get(i).norm;
			}		
			
		}
		return unmat;
	}
	
	public void setPRMScores(int charge){
		for( int i=0; i<this.size(); i++){
			setPRMScore( i, charge, 0.9 );
		}
	}
	
	public PNode getMatchedPNode( double mz ){
		double it=0;
		PNode anno = null;
		int id= binarySearch( mz- Mutables.fragmentTolerance );
		if( this.get(id).getMass() < (mz - Mutables.fragmentTolerance) ) return anno;
		
		while( this.get(id).getMass() <= mz + Mutables.fragmentTolerance ){
			if( this.get(id).getNorm() > it ){
				it = this.get(id).getNorm();
				anno = this.get(id);
			}
			id++;
			if( id == this.size() )
				break;
		}
		return anno;
	}
	
	public double getNormOfMatchedPeak( double mz ){
		int index;
		if( (index= getPosition(mz)) == -1 ) return 0;
		else return this.get(index).getNorm();
	}
	public double getNormOfMatchedPeak( double mz, int rank ){
		double it=0;
		int anno = -1;
		int id= binarySearch( mz- Mutables.fragmentTolerance );
		if( this.get(id).getMass() < (mz - Mutables.fragmentTolerance) ) return it;

		while( this.get(id).getMass() <= mz + Mutables.fragmentTolerance ){
			if( this.get(id).getNorm() > it ){
				it = this.get(id).getNorm();
				anno = id;
			}
			id++;
		
			if( id == this.size() )
				break;
		}
		if( anno != -1 && !this.get(anno).isAssigned() && this.get(anno).localRank <= rank ) {
			this.get(anno).assign(true);
			return it;
		}
		else return 0;

	}
	
	public int getIndexOfMatchedPeak( double mz ){
		return getPosition(mz);
	}
	
	public int getIndexOfMatchedPeak(double mz, String type, double peakWidth){
		double it=0;
		int anno = -1;
		int id= binarySearch( mz- peakWidth );
		if( this.get(id).getMass() < (mz - peakWidth) ) return anno;
		
		while( this.get(id).getMass() <= mz + peakWidth ){
			if( this.get(id).getIntensity() > it ){
				it = this.get(id).getIntensity();
				anno = id;
			}
			id++;
			if( id == this.size() )
				break;
		}		
		if( anno != -1 && !this.get(anno).isAssigned() ){
			this.get(anno).assign(true, type);
		}
		return anno;
	}
	public int getIndexOfMatchedPeak(double mz, double err, String type){
		double it=0;
		int anno = -1;
		int id= binarySearch( mz- Mutables.fragmentTolerance );
		if( this.get(id).getMass() < (mz - Mutables.fragmentTolerance) ) return anno;
		
		while( this.get(id).getMass() <= mz+err+ Mutables.fragmentTolerance ){
			if( this.get(id).getIntensity() > it ){
				it = this.get(id).getIntensity();
				anno = id;
			}
			id++;
			if( id == this.size() )
				break;
		}		
		if( anno != -1 && !this.get(anno).isAssigned() ){
			this.get(anno).assign(true, type);
		}
		return anno;
	}
	
	private int binarySearch( double left )
	{
		int index;	
		if( left <= this.get(0).getMass() )
			index= 0;
		else if( left > this.get(this.size()-1).getMass() )
			index= this.size()-1;
		else
		{
			int M, L= 0, R= this.size()-1;
			while( R - L > 1 )
			{
				M= ( L + R ) /2;

				if( left <= this.get(M).getMass() )
					R= M;
				else
					L= M;
			}
			index= R;
		}	
		return index;
	}
	
	private int getPosition( double mz ){
		double delta=Mutables.fragmentTolerance+1;
		int anno = -1;
		int id= binarySearch( mz- Mutables.fragmentTolerance );
		if( this.get(id).getMass() < (mz - Mutables.fragmentTolerance) ) {
			return anno;
		}
		while( this.get(id).getMass() <= mz + Mutables.fragmentTolerance ){
			if( Math.abs(this.get(id).mass-mz) < delta ){
				delta = Math.abs(this.get(id).mass - mz);
				anno = id;
			}
			id++;
			if( id == this.size() )
				break;
		}	
		
		if( anno == -1 || this.get(anno).isAssigned() ) {
			return -1;
		}
		else {
			this.get(anno).assign(true);
			return anno;
		}
	}

	public String toString(){
		StringBuffer a = new StringBuffer();
		for(PNode p : this){
		
			if( p.assigned ){
				if( p.annotation.startsWith("DB2") ){
					if( p.localRank < 3 ) a.append(p +"\r\n");
				}
				else if( p.annotation.startsWith("CS") || p.annotation.startsWith("SS") ){
					if( p.localRank < 5 ) a.append(p +"\r\n");
				}
				else if( p.localRank < 6 )
					a.append(p +"\r\n");
			}
		}
		return a.toString();
	}
		
	public double getSupportingScore( int peak, int type, int charge, double supportingFactor ){
		
		double isotopeDelta = Constants.IsotopeSpace/charge, NH3Delta = Constants.NH3/charge, H2ODelta = Constants.H2O/charge;
		
		double score= 0, targetmz;
		double ionMZ = this.get(peak).getMass(), ionIT = this.get(peak).getIntensity();			
		
		// check this is isotope???
	//	int OK = -1;    	
		double okmax = 0;
		targetmz= ionMZ-isotopeDelta;
		for(int i=peak-1; i>-1; i--){
			if( this.get(i).getMass() < targetmz-Mutables.massToleranceForDenovo ) break;			
			else if( Math.abs(this.get(i).getMass()-targetmz) < Mutables.massToleranceForDenovo ){
				if( this.get(i).getIntensity() > okmax ){
					okmax = this.get(i).getIntensity();
				//	OK = i;
				}
			}
		}
		if( okmax > ionIT/isotopeRatio ) return 0;
		
		//start finding supporting peaks
		int ISO = peak + 1; // plus isotope peak;    
		double prevISO = ionIT;
		boolean isoDecent = false;
		targetmz= ionMZ + isotopeDelta;
		for( int nth_iso=1 ;  ; nth_iso++  ){
			int H = -1;
			double imax = 0;				
			for(int i=ISO; i<this.size(); i++){
				if( this.get(i).getMass() > targetmz+Mutables.massToleranceForDenovo ) break;			
				else if( Math.abs(this.get(i).getMass()-targetmz) < Mutables.massToleranceForDenovo ){
					if( this.get(i).getIntensity() > imax ) {
						imax = this.get(i).getIntensity();
						H = i;
					}
				}
			}
			if( H == -1 || this.get(H).isAssigned() ) break;
			if( prevISO < imax/isotopeRatio ) break;
			if( isoDecent && prevISO < imax ) break;
			if( prevISO > imax ) isoDecent = true;
				
			score += this.get(H).getNorm()/nth_iso;	
			this.get(H).assign(true);
			ISO = H+1;
			targetmz= this.get(H).getMass() + isotopeDelta;		
			prevISO= imax;
		}//*/

		int NLOSS = -1;    	
		double lossmax=0, lossmz=0;
		targetmz= ionMZ-H2ODelta;
		for(int i=peak-1; i>-1; i--){
			if( this.get(i).getMass() < targetmz-Mutables.massToleranceForDenovo ) break;			
			else if( Math.abs(this.get(i).getMass()-(ionMZ-NH3Delta)) < Mutables.massToleranceForDenovo ){
				if( this.get(i).getIntensity() > lossmax ) {
					lossmax = this.get(i).getIntensity();
					lossmz = this.get(i).getMass();
					NLOSS = i;
				}
			}
			else if( Math.abs(this.get(i).getMass()-targetmz) < Mutables.massToleranceForDenovo ){
				if( this.get(i).getIntensity() > lossmax ) {
					lossmax = this.get(i).getIntensity();
					lossmz = this.get(i).getMass();
					NLOSS = i;
				}
			}
		}
		
		if( NLOSS != -1 && !this.get(NLOSS).isAssigned() ){
			if( type == 0 || lossmax < ionIT )
			{
				score += this.get(NLOSS).getNorm()*supportingFactor;
				this.get(NLOSS).assign(true);
				
				int NL_H = -1; // plus isotope peak    	
				double nsmax = 0;
				targetmz= lossmz+isotopeDelta;
				for(int i=NLOSS+1; i<this.size(); i++){
					if( this.get(i).getMass() > targetmz+Mutables.massToleranceForDenovo ) break;			
					else if( Math.abs(this.get(i).getMass()-targetmz) < Mutables.massToleranceForDenovo ){
						if( this.get(i).getIntensity() > nsmax ) {
							nsmax = this.get(i).getIntensity();
							NL_H = i;
						}
					}
				}
				if( NL_H != -1 && !this.get(NL_H).isAssigned() ) {
					score += this.get(NL_H).getNorm()*supportingFactor;	
					this.get(NL_H).assign(true);
				}
			}			
		}//*/
		
		return score;
	}
	
	public void setPRMScore( int peak, int charge, double supportingFactor ){
		
		this.get(peak).b_prm = this.get(peak).y_prm = 0;
		double isotopeDelta = Constants.IsotopeSpace/charge, NH3Delta = Constants.NH3/charge, H2ODelta = Constants.H2O/charge;
	
		double targetmz;
		double score = this.get(peak).getNorm(), ionMZ= this.get(peak).getMass(), ionIT= this.get(peak).getIntensity();			
		
		// check this is isotope???
		int OK = -1;    	
		double okmax = 0;
		targetmz= ionMZ-isotopeDelta;
		for(int i=peak-1; i>-1; i--){
			if( this.get(i).getMass() < targetmz-Mutables.massToleranceForDenovo ) break;			
			else if( Math.abs(this.get(i).getMass()-targetmz) < Mutables.massToleranceForDenovo ){
				if( this.get(i).getIntensity() > okmax ) {
					okmax = this.get(i).getIntensity();
					OK = i;
				}
			}
		}
		if( OK != -1 ) {
			if( okmax > ionIT ) this.get(peak).b_prm = this.get(peak).y_prm = 0;
			else this.get(peak).b_prm = this.get(peak).y_prm = score;
			return;
		}//*/
		
		int ISO = peak + 1; // plus isotope peak    
		double prevISO = ionIT;
		boolean isoDecent = false;
		targetmz= ionMZ + isotopeDelta;
		for( int nth_iso=1 ;  ; nth_iso++  ){
			int H = -1;
			double imax = 0;				
			for(int i=ISO; i<this.size(); i++){
				if( this.get(i).getMass() > targetmz+Mutables.massToleranceForDenovo ) break;			
				else if( Math.abs(this.get(i).getMass()-targetmz) < Mutables.massToleranceForDenovo ){
					if( this.get(i).getIntensity() > imax ) {
						imax = this.get(i).getIntensity();
						H = i;
					}
				}
			}
		//	if( H == -1 ) break;
			if( H == -1 || prevISO < imax/isotopeRatio ) break;
			if( isoDecent && prevISO < imax ) break;
			if( prevISO > imax ) isoDecent = true;
			
			score += this.get(H).getNorm()/nth_iso;	
			ISO = H+1;
			targetmz= this.get(H).getMass() + isotopeDelta;		
			prevISO = imax;
		}//*/
		
		this.get(peak).b_prm = this.get(peak).y_prm = score;

		int NLOSS = -1;    	
		double lossmax=0, lossmz=0;
		targetmz= ionMZ-H2ODelta;
		for(int i=peak-1; i>-1; i--){
			if( this.get(i).getMass() < targetmz-Mutables.massToleranceForDenovo ) break;			
			else if( Math.abs(this.get(i).getMass()-(ionMZ-NH3Delta)) < Mutables.massToleranceForDenovo ){
				if( this.get(i).getIntensity() > lossmax ) {
					lossmax = this.get(i).getIntensity();
					lossmz = this.get(i).getMass();
					NLOSS = i;
				}
			}
			else if( Math.abs(this.get(i).getMass()-targetmz) < Mutables.massToleranceForDenovo ){
				if( this.get(i).getIntensity() > lossmax ) {
					lossmax = this.get(i).getIntensity();
					lossmz = this.get(i).getMass();
					NLOSS = i;
				}
			}
		}
		
		if( NLOSS != -1 ){
			double lossScore = this.get(NLOSS).getNorm();		
			int NL_H = -1; // plus isotope peak    	
			double nsmax = 0;
			targetmz= lossmz+isotopeDelta;
			for(int i=NLOSS+1; i<this.size(); i++){
				if( this.get(i).getMass() > targetmz+Mutables.massToleranceForDenovo ) break;			
				else if( Math.abs(this.get(i).getMass()-targetmz) < Mutables.massToleranceForDenovo ){
					if( this.get(i).getIntensity() > nsmax ) {
						nsmax = this.get(i).getIntensity();
						NL_H = i;
					}
				}
			}
			if( NL_H != -1 ) lossScore += this.get(NL_H).getNorm();
			
			this.get(peak).b_prm += lossScore*supportingFactor;
			if( lossmax < ionIT ) this.get(peak).y_prm += lossScore*supportingFactor;		
		}//*/
	}

	public void sortRevByPNodeIntensity(){
		this.sort(Collections.reverseOrder(new PNodeIntComparator()));
	}
	

	public double getSupportingScore( int peak, int type ){
		double score = 0;	
		double ionMZ= this.get(peak).getMass(), ionIT= this.get(peak).getIntensity();			
		
		// check this is isotope???
		int OK = -1;    	
		double okmax = 0;
		for(int i=peak-1; i>-1; i--){
			if( this.get(i).getMass() < ionMZ-1-Mutables.massToleranceForDenovo ) break;			
			else if( Math.abs(this.get(i).getMass()- (ionMZ-1)) < Mutables.massToleranceForDenovo ){
				if( this.get(i).getIntensity() > okmax ) {
					okmax = this.get(i).getIntensity();
					OK = i;
				}
			}
		}
		if( OK != -1 ) return 0;
				
		int H = -1; // plus isotope peak    	
		double imax = 0;
		for(int i=peak+1; i<this.size(); i++){
			if( this.get(i).getMass() > ionMZ+1+Mutables.massToleranceForDenovo ) break;			
			else if( Math.abs(this.get(i).getMass()- (ionMZ+1)) < Mutables.massToleranceForDenovo ){
				if( this.get(i).getIntensity() > imax ) {
					imax = this.get(i).getIntensity();
					H = i;
				}
			}
		}
		if( H != -1 && !this.get(H).isAssigned() ) {
			score += this.get(H).getNorm();	
			this.get(H).assign(true);
		}
		
		int H2O = -1, NH3 = -1;    	
		double hmax=0, nmax=0;
		for(int i=peak-1; i>-1; i--){
			if( this.get(i).getMass() < ionMZ-Constants.H2O-Mutables.massToleranceForDenovo ) break;			
			else if( Math.abs(this.get(i).getMass()- (ionMZ-Constants.NH3)) < Mutables.massToleranceForDenovo ){
				if( this.get(i).getIntensity() > nmax ) {
					nmax = this.get(i).getIntensity();
					NH3 = i;
				}
			}
			else if( Math.abs(this.get(i).getMass()- (ionMZ-Constants.H2O)) < Mutables.massToleranceForDenovo ){
				if( this.get(i).getIntensity() > hmax ) {
					hmax = this.get(i).getIntensity();
					H2O = i;
				}
			}
		}	
		
		if( NH3 != -1 && !this.get(NH3).isAssigned() ){
			if( type == 1 ){
				score += this.get(NH3).getNorm();
				this.get(NH3).assign(true);
			}
			else if( nmax < ionIT ){
				score += this.get(NH3).getNorm();	
				this.get(NH3).assign(true);
			}
		}
		if( H2O != -1 && !this.get(H2O).isAssigned() ){
			if( type == 1 ){
				score += this.get(H2O).getNorm();
				this.get(H2O).assign(true);
			}
			else if( hmax < ionIT ){
				score += this.get(H2O).getNorm();
				this.get(H2O).assign(true);
			}
		}//
		return score;
	}
}







