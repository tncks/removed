package modi;

import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.Collections;
import java.util.StringTokenizer;

public class ScanCap implements Comparable<ScanCap> {
	private String 	title;
	private double 	pmz;
	private double 	neutralMW;
	private int 	charge;
	private int 	scanNo;
	private long 	offset;
	
	private static double tolerance= Constants.massToleranceForDenovo;
	
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
	
	public Spectrum getSpectrum( RandomAccessFile in ) throws IOException {
		
		if( Constants.rangeForIsotopeIncrement != 0 ){
			Constants.maxNoOfC13 = (int)Math.ceil( neutralMW / Constants.rangeForIsotopeIncrement );			
		}
		
		if( Constants.PPMTolerance != 0 ) {
			Constants.precursorAccuracy = Constants.PPMtoDalton(neutralMW, Constants.PPMTolerance);
		}
		Constants.precursorTolerance = Constants.precursorAccuracy + Constants.maxNoOfC13*Constants.IsotopeSpace;
	
		String s;
		
		ArrayList<RawP> rawPL = new ArrayList<RawP>();
		in.seek( this.offset );
		while( (s = in.readLine()) != null ) {
			StringTokenizer token = new StringTokenizer(s);
			if( token.countTokens() > 1 ){
				if( !Character.isDigit(s.charAt(0)) ) break;
				rawPL.add( new RawP(Double.parseDouble(token.nextToken()), Double.parseDouble(token.nextToken())) );
			}
			else break;
		}
		Collections.sort( rawPL );
		
	//	processingiTRAQ(rawPL);
		
		int index = 0;
		Spectrum spectrum = new Spectrum( this.pmz, this.charge, this.title );
		
		double basePeakIntensity=0, TIC=0;
		double tarMass=0, tarInten=0;
		for( RawP rp : rawPL ) {
			double mass = rp.mz;
			double intensity = rp.it;
			
			if( intensity <= 0 || mass <= 0 ) continue;
			if( mass > neutralMW ) continue;
		//	if( Math.abs( mass-pmz ) < 2. ) continue;	
			
			if( ( mass - tarMass ) < tolerance ){
				double sum = tarInten + intensity;
				tarMass = tarMass*(tarInten/sum)+ mass*(intensity/sum);
				tarInten += intensity;
				spectrum.get(index-1).set(tarMass, tarInten);
			}
			else{
				spectrum.add( new Peak(index++, mass, intensity) );
				tarMass = mass;
				tarInten = intensity; 
			}
			TIC += intensity;
			if( tarInten > basePeakIntensity )
				basePeakIntensity= tarInten;
		}		
		spectrum.setExtraInformation( basePeakIntensity, TIC );
		
		Constants.gapTolerance = Constants.fragmentTolerance*2;
		Constants.nonModifiedDelta = (Constants.precursorTolerance<Constants.massToleranceForDenovo)? Constants.precursorTolerance : Constants.massToleranceForDenovo;
				
		if( Constants.precursorTolerance > Constants.gapTolerance )
			Constants.gapTolerance = Constants.precursorTolerance;		
		
		return spectrum; 
	}
	
	private class RawP implements Comparable<RawP> {
		double mz;
		double it;
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
	
	private void xprocessingiTRAQ(ArrayList<RawP> rawPL){
		
		ArrayList<RawP> reporters = new ArrayList<RawP>();
		reporters.add( new RawP(114.1105, Constants.fragmentTolerance) );//reporterIon114
		reporters.add( new RawP(115.1074, Constants.fragmentTolerance) );//reporterIon115
		reporters.add( new RawP(116.1107, Constants.fragmentTolerance) );//reporterIon116
		reporters.add( new RawP(117.1141, Constants.fragmentTolerance) );//reporterIon117
		
		reporters.add( new RawP(Constants.NTERM_FIX_MOD + Constants.Proton, Constants.fragmentTolerance) );//iTRAQ TAG
		
		int fragCS = 1;
		while( true ){
			double compItraqTag = (this.neutralMW - Constants.NTERM_FIX_MOD + Constants.Proton*fragCS)/fragCS;
	
			double secondIso = compItraqTag + Constants.IsotopeSpace/fragCS;
			double thirdIso  = secondIso + Constants.IsotopeSpace/fragCS;
			double forthIso  = thirdIso + Constants.IsotopeSpace/fragCS;
			
			reporters.add( new RawP(compItraqTag, Constants.fragmentTolerance) );//precursor without iTRAQ TAG
			reporters.add( new RawP(secondIso, Constants.fragmentTolerance) );//precursor without iTRAQ TAG
			reporters.add( new RawP(thirdIso, Constants.fragmentTolerance) );//precursor without iTRAQ TAG
			reporters.add( new RawP(forthIso, Constants.fragmentTolerance) );//precursor without iTRAQ TAG//*/

//			reporters.add( new RawP(compItraqTag, (3*Constants.IsotopeSpace/fragCS)+Constants.fragmentTolerance) );//precursor without iTRAQ TAG
			for(int i=1; i<=Constants.maxNoOfC13; i++){
				reporters.add( new RawP(compItraqTag-i*Constants.IsotopeSpace/fragCS, Constants.fragmentTolerance) );//precursor without iTRAQ TAG
			}//*/

			if( ++fragCS >= this.charge ) break;
		}
		
		Collections.sort(reporters);
		
		int start = 0;
		for( RawP rp : reporters ){
			for (int i=start; i<rawPL.size(); i++){
				if( rawPL.get(i).mz < rp.mz-rp.it ) continue;
				else if( rawPL.get(i).mz > rp.mz+rp.it ) {
					start = i;
					break;
				}				
				rawPL.remove(i);
				i--;
			}
		}
	
	}
}
