package msutil;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;

import modi.AnsPeptide;
import modi.Constants;
import modi.PTM;
import modi.PTMDB;
import modi.ScanCap;

import org.jdom.Attribute;
import org.jdom.Element;

import processedDB.PeptideMatchToProtein;
import processedDB.ProtDatabase;
import scaniter.MSMScan;

public class Modxml {
	public static void startWriting( PrintStream xml ){
		xml.println( "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" );
		xml.println("<msms_search>");
	}
	public static void endWriting( PrintStream xml ){
		xml.println( "</msms_search>" );
	}
	public static Element get_search_summary(){
		
		Element summary = new Element( "summary" );
		summary.setAttribute( new Attribute( "engine", Constants.engine ) );
		summary.setAttribute( new Attribute( "version", Constants.engineVersion ) );
		summary.setAttribute( new Attribute( "date", Constants.runDate ) );
		summary.setAttribute( new Attribute( "user", Constants.runUser ) );
		summary.setAttribute( new Attribute( "title", Constants.runTitle ) );
		
		Element dataset = new Element( "dataset" );
		dataset.setAttribute( new Attribute( "path", Constants.SPECTRUM_LOCAL_PATH ) );
		dataset.setAttribute( new Attribute( "format", String.valueOf( Constants.SPECTRA_FILE_TYPE ) ) );
		dataset.setAttribute( new Attribute( "instrument", Constants.INSTRUMENT_NAME ) );
		summary.addContent( dataset );
		
		Element database = new Element( "database" );
		database.setAttribute(new Attribute( "path", Constants.PROTEIN_DB_LOCAL_PATH ));
		summary.addContent( database );
		
		Element enzyme = new Element( "enzyme" );
		enzyme.setAttribute(new Attribute( "name", Constants.protease.getName() ));
		enzyme.setAttribute(new Attribute( "cleave", Constants.protease.getCleave() ));
		enzyme.setAttribute(new Attribute( "term", Constants.protease.getTerm() ));
		summary.addContent( enzyme );
		
		Element parameters = new Element( "parameters" );
		Element enzyme_constraint = new Element( "enzyme_constraint" );
		enzyme_constraint.setAttribute( new Attribute( "max_miss_cleavages", String.valueOf(Constants.missCleavages) ) );
		enzyme_constraint.setAttribute( new Attribute( "min_number_termini", String.valueOf(Constants.numberOfEnzymaticTermini) ) );
		parameters.addContent( enzyme_constraint );
		
		Element peptide_mass_tol = new Element( "peptide_mass_tol" );
		if( Constants.PPMTolerance == 0 ){
			peptide_mass_tol.setAttribute( new Attribute( "value", String.format("%.4f", Constants.precursorTolerance) ) );
			peptide_mass_tol.setAttribute( new Attribute( "unit", "Da" ) );	
		}
		else{
			peptide_mass_tol.setAttribute( new Attribute( "value", String.format("%.1f", Constants.PPMTolerance) ) );
			peptide_mass_tol.setAttribute( new Attribute( "unit", "ppm" ) );	
		}
			
		parameters.addContent( peptide_mass_tol );
		
		Element fragment_ion_tol = new Element( "fragment_ion_tol" );
		fragment_ion_tol.setAttribute( new Attribute( "value", String.format("%.4f", Constants.fragmentTolerance) ) );
		fragment_ion_tol.setAttribute( new Attribute( "unit", "Da" ) );
		parameters.addContent( fragment_ion_tol );
		
		Element modified_mass_range = new Element( "modified_mass_range" );	
		modified_mass_range.setAttribute( new Attribute( "min_value", String.format("%.1f", Constants.minModifiedMass) ) );
		modified_mass_range.setAttribute( new Attribute( "max_value", String.format("%.1f", Constants.maxModifiedMass) ) );
		parameters.addContent( modified_mass_range );

		summary.addContent( parameters );
		
		Element modifications = new Element( "modifications" );
		Element variable = new Element( "variable" );
		if( Constants.variableModifications.size() > 0 ){
			variable.setAttribute( new Attribute( "num", String.valueOf(Constants.variableModifications.size()) ) );
			for( PTM p : Constants.variableModifications ){
				Element mod = new Element( "mod" );
				mod.setAttribute( new Attribute( "name", p.getName() ) );
				mod.setAttribute( new Attribute( "site", p.getSite() ) );
				mod.setAttribute( new Attribute( "position", String.valueOf(p.getPTMPosition()) ) );
				mod.setAttribute( new Attribute( "massdiff", String.format("%.4f", p.getMassDifference()) ) );
				variable.addContent( mod );
			}
		}
		modifications.addContent( variable );
		Element fixed = new Element( "fixed" );
		if( Constants.fixedModifications.size() > 0 ){
			fixed.setAttribute( new Attribute( "num", String.valueOf(Constants.fixedModifications.size()) ) );
			for( PTM p : Constants.fixedModifications ){
				Element mod = new Element( "mod" );
				mod.setAttribute( new Attribute( "name", p.getName() ) );
				if( p.getSite().compareToIgnoreCase("N-TERM") == 0 ) mod.setAttribute( new Attribute( "site", "[" ) );
				else if( p.getSite().compareToIgnoreCase("C-TERM") == 0 ) mod.setAttribute( new Attribute( "site", "]" ) );
				else mod.setAttribute( new Attribute( "site", p.getSite() ) );
				mod.setAttribute( new Attribute( "position", String.valueOf(p.getPTMPosition()) ) );
				mod.setAttribute( new Attribute( "massdiff", String.format("%.4f", p.getMassDifference()) ) );
				fixed.addContent( mod );
			}
		}
		modifications.addContent( fixed );
		
		summary.addContent( modifications );	
		return summary;	
	}
	
	public static Element get_spectrum_query_element( int index, ScanCap spectrum, ProtDatabase protDB, PTMDB ptmDB, ArrayList<AnsPeptide> candidates ){
		//org modi
		Element spectrum_guery = new Element( "query" );
		spectrum_guery.setAttribute(new Attribute( "id", String.valueOf(index) ));
		spectrum_guery.setAttribute(new Attribute( "spectrum", spectrum.getTitle() ));
		spectrum_guery.setAttribute(new Attribute( "observed_mz", String.format("%.4f", spectrum.getPMZ()) ));
		spectrum_guery.setAttribute(new Attribute( "cs", String.valueOf(spectrum.getCharge()) ));
		spectrum_guery.setAttribute(new Attribute( "position", String.valueOf(spectrum.getOffset()) ));
		
		for( int k=0; k<candidates.size(); k++ ){					
			spectrum_guery.addContent( candidates.get(k).get_peptide_hit_element(k+1, protDB, spectrum.getObservedMW(), ptmDB) );
		}			
		return spectrum_guery;		
	}
	
	public static Element get_spectrum_query_element( MSMScan spectrum, PTMDB ptmDB, ArrayList<AnsPeptide> candidates, 
			HashMap<String, ArrayList<PeptideMatchToProtein>> seqToProtMap ){
		
		Element spectrum_guery = new Element( "query" );
		spectrum_guery.setAttribute(new Attribute( "id", String.valueOf(spectrum.getSpecIndex()) ));
	//	spectrum_guery.setAttribute(new Attribute( "spectrum", String.valueOf(spectrum.getSpecIndex()) ));
		spectrum_guery.setAttribute(new Attribute( "spectrum", spectrum.getTitle() ));
		spectrum_guery.setAttribute(new Attribute( "observed_mz", String.format("%.4f", spectrum.getPMZ()) ));
		spectrum_guery.setAttribute(new Attribute( "cs", String.valueOf(spectrum.getCharge()) ));
		spectrum_guery.setAttribute(new Attribute( "position", String.valueOf(spectrum.getOffset()) ));
		
		for( int k=0; k<candidates.size(); k++ ){					
			spectrum_guery.addContent( candidates.get(k).get_peptide_hit_element(k+1, spectrum.getObservedMW(), ptmDB, seqToProtMap.get(candidates.get(k).getPeptideSequence())) );
		}			
		return spectrum_guery;		
	}
	
}
