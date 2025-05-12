package scaniter;

import java.util.ArrayList;
import java.util.Collections;
import moda.MODaConst;

import modi.Constants;
import modi.Mutables;
import modi.Peak;
import modi.Spectrum;
import msutil.MSMass;

public class MSMScan {
	

	static final double minMW = 8 * MSMass.getMinAAMass() + 18;
	
	private final String 		title;
	private final int 		specIndex;
	private final int 		scanNo;
	private final double 		pmz;
	private final double 		neutralMW;
	private final int 		charge;
	private Spectrum 	peaklist;

	public double getPrecursorTolerance() {
		return precursorTolerance;
	}

	public double getPrecursorAccuracy() {
		return precursorAccuracy;
	}

	public double getGapTolerance() {
		return gapTolerance;
	}

	public double getNonModifiedDelta() {
		return nonModifiedDelta;
	}

	public int getMaxNoOfC13() {
		return maxNoOfC13;
	}

	private double 	precursorTolerance = 0;
	private double 	precursorAccuracy= 0;
	private double 	gapTolerance = 0;
	private double 	nonModifiedDelta = 0;
	private int		maxNoOfC13 = 0;


	public MSMScan(int index, double pmz, int charge){
		this.title 		= "";
		this.specIndex	= index;
		this.scanNo		= 0;	
		this.pmz 		= pmz;
		this.charge 	= charge;
		this.neutralMW 	= (pmz - Constants.Proton)*charge;
	}
	
	public MSMScan(String title, int index, int sn, double pmz, int charge){		
		this.title 		= title;
		this.specIndex  = index;
		this.scanNo		= sn;
		this.pmz 		= pmz;
		this.charge 	= charge;
		this.neutralMW 	= (pmz - Constants.Proton)*charge;
	}

	public Spectrum getSpectrum() {
		return peaklist; 
	}

	public boolean setSpectrum(ArrayList<RawPeak> rawPL) {
		
		if( neutralMW < minMW || Constants.maxPeptideMass < neutralMW ) return false;
				
		if( Mutables.reporterMassOfIsobaricTag != null ) removeReporterIons(rawPL, Mutables.reporterMassOfIsobaricTag);
		
		if( Constants.rangeForIsotopeIncrement != 0 ) maxNoOfC13 = (int)Math.ceil( neutralMW / Constants.rangeForIsotopeIncrement );
		else maxNoOfC13 = Mutables.DEFAULT_MAXNOOFC13; ///////////////////// default ?
		
		if( Constants.PPMTolerance != 0 ) precursorAccuracy = Constants.PPMtoDalton( neutralMW, Constants.PPMTolerance );
		else precursorAccuracy = Mutables.DEFAULT_PRECURSORACCURACY; ///////////////////// default ?
		
		precursorTolerance = precursorAccuracy + maxNoOfC13*Constants.IsotopeSpace;
		
		int index = 0;
		Spectrum spectrum = new Spectrum( this.pmz, this.charge, this.title );
		
		double basePeakIntensity=0, TIC=0;
		double tarMass=0, tarInten=0;
		for( RawPeak rp : rawPL ) {		
			double mass = rp.mz;
			double intensity = rp.it;
			if( intensity <= 0 || mass <= 0 ) continue;
			if( mass > neutralMW ) continue;
			
			if( ( mass - tarMass ) < Mutables.massToleranceForDenovo ){
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
		
		gapTolerance = Mutables.fragmentTolerance*2;
		nonModifiedDelta = (precursorTolerance < Mutables.massToleranceForDenovo)? precursorTolerance : Mutables.massToleranceForDenovo;
				
		if( precursorTolerance > gapTolerance ) gapTolerance += precursorTolerance;
		
		if( spectrum.size() < MODaConst.minPeaksCount ) peaklist = null;
		else peaklist = spectrum;
		
		return (peaklist!=null);
	}
	
	private void removeReporterIons( ArrayList<RawPeak> rawPL, double[] removedMasses ){
	
		ArrayList<RawPeak> reporters = new ArrayList<>();
		for(int i=1; i<removedMasses.length; i++)
			reporters.add( new RawPeak(removedMasses[i], Mutables.fragmentTolerance) );
		
		reporters.add( new RawPeak(removedMasses[0] + Constants.Proton, Mutables.fragmentTolerance) );

		int fragCS = 1;
		while( true ){
			double compItraqTag = (this.neutralMW - removedMasses[0] + Constants.Proton*fragCS)/fragCS;
	
			double secondIso = compItraqTag + Constants.IsotopeSpace/fragCS;
			double thirdIso  = secondIso + Constants.IsotopeSpace/fragCS;
			double forthIso  = thirdIso + Constants.IsotopeSpace/fragCS;
			
			reporters.add( new RawPeak(compItraqTag, Mutables.fragmentTolerance) );
			reporters.add( new RawPeak(secondIso, Mutables.fragmentTolerance) );
			reporters.add( new RawPeak(thirdIso, Mutables.fragmentTolerance) );
			reporters.add( new RawPeak(forthIso, Mutables.fragmentTolerance) );

			for(int i=1; i<=maxNoOfC13; i++){
				reporters.add( new RawPeak(compItraqTag-i*Constants.IsotopeSpace/fragCS, Mutables.fragmentTolerance) );
			}
			if( ++fragCS >= this.charge ) break;
		}

		Collections.sort(reporters);
		
		int start = 0;
		for( RawPeak rp : reporters ){
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

	public double 	getObservedMW(){ return neutralMW; }

	public String 	getHeader(){ return String.format("%d\t%.4f\t%d\t%d\t%s",
			specIndex, neutralMW, charge, scanNo, title); }
	
}

























