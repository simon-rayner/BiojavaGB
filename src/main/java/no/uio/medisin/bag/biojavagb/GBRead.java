/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.biojavagb;


import java.io.File;
import java.io.FileInputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;


import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava.nbio.core.sequence.compound.DNACompoundSet;
 
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.sequence.features.FeatureInterface;
import org.biojava.nbio.core.sequence.features.Qualifier;
import org.biojava.nbio.core.sequence.io.DNASequenceCreator;
import org.biojava.nbio.core.sequence.io.GenbankReader;
import org.biojava.nbio.core.sequence.io.GenbankReaderHelper;
import org.biojava.nbio.core.sequence.io.GenericGenbankHeaderParser;
import org.biojava.nbio.core.sequence.io.ProteinSequenceCreator;
import org.biojava.nbio.core.sequence.loader.GenbankProxySequenceReader;
import org.biojava.nbio.core.sequence.template.AbstractSequence;

import org.apache.commons.cli.Options;
import org.apache.commons.lang3.StringUtils;
//import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;




/**
 * linux
 * @author simonray
 */

public class GBRead {
 
    static Options options = new Options();
    static Logger logger = LogManager.getRootLogger();
    
    public static void main(String[] args) throws Exception{
        
        logger.info("\n\n" + StringUtils.repeat("*", 80) + "\n");
        logger.info("                  GenBank parse");
        logger.info("\n\n" + StringUtils.repeat("*", 80) + "\n");
//        String ax = StringUtils.repeat("*", 80);
//        String gx = Strings.repeat("*", 3);
        
	/*
	 * Method 1: With the GenbankProxySequenceReader
	 */
	//Try with the GenbankProxySequenceReader
	GenbankProxySequenceReader<AminoAcidCompound> genbankProteinReader 
	= new GenbankProxySequenceReader<AminoAcidCompound>("/tmp", "NP_000257", AminoAcidCompoundSet.getAminoAcidCompoundSet());
	ProteinSequence proteinSequence = new ProteinSequence(genbankProteinReader);
	genbankProteinReader.getHeaderParser().parseHeader(genbankProteinReader.getHeader(), proteinSequence);
	logger.info("Sequence" + "(" + proteinSequence.getAccession() + "," + proteinSequence.getLength() + ")=" +
            proteinSequence.getSequenceAsString().substring(0, 10) + "...\n\n");
 
	GenbankProxySequenceReader<NucleotideCompound> genbankDNAReader 
            = new GenbankProxySequenceReader<NucleotideCompound>("/tmp", "NM_001126", DNACompoundSet.getDNACompoundSet());
	DNASequence dnaSequence = new DNASequence(genbankDNAReader);
	genbankDNAReader.getHeaderParser().parseHeader(genbankDNAReader.getHeader(), dnaSequence);
	logger.info("Sequence" + "(" + dnaSequence.getAccession() + "," + dnaSequence.getLength() + ")=" +
            dnaSequence.getSequenceAsString().substring(0, 10) + "...\n\n");
	/*
	 * Method 2: With the GenbankReaderHelper
	 */
	//Try with the GenbankReaderHelper
	File dnaFile = new File("src/test/resources/NM_000266.gb");		
	File protFile = new File("src/test/resources/BondFeature.gb");
 
        logger.info("\n\nREAD dna file " + dnaFile + "\n\n");
	LinkedHashMap<String, DNASequence> dnaSequences = GenbankReaderHelper.readGenbankDNASequence( dnaFile );
	for (DNASequence sequence : dnaSequences.values()) {
	    	//logger.info( sequence.getSequenceAsString() );

                String notes = sequence.getOriginalHeader();
                sequence.getDescription();
                List fr = sequence.getFeaturesByType("source");
                
//                List<FeatureInterface<AbstractSequence<NucleotideCompound>, NucleotideCompound>> fr = sequence.getFeaturesByType("SOURCE");
//                for (FeatureInterface fi : fr) {
//                    logger.info("DESCRIPTION\t" + fi.getDescription());                    
//                }                
                logger.info(sequence.getDatabaseReferences());
                logger.info(sequence.getDescription());

                List<FeatureInterface<AbstractSequence<NucleotideCompound>, NucleotideCompound>> fl =   sequence.getFeatures();
                for (FeatureInterface fi : fl) {
                    logger.info("DESCRIPTION\t" + fi.getDescription());
                    logger.info("CHILDREN\t" + fi.getChildrenFeatures());
                    logger.info("LOCATIONS\t" + fi.getLocations());
                    logger.info("PARENT\t" + fi.getParentFeature());
                
                    logger.info("QUALIFIERS");
                    HashMap <String, Qualifier> quals = fi.getQualifiers();
                    for(Map.Entry<String, Qualifier> entry : quals.entrySet()){
                        logger.info("--\t" + entry.getKey() + "\t|\t" + entry.getValue().getName() 
                                + "  /  " + entry.getValue().getValue() + "\\" + entry.getValue().toString());                       
                    }

                    logger.info("\n");
                    logger.info("SHORT\t" + fi.getShortDescription());
                    logger.info("SOURCE\t" + fi.getSource());
                    logger.info("TYPE\t" + fi.getType());
                    logger.info("HASHCODE\t" + fi.hashCode());
                    logger.info("-");
                }
//                Annotation seqAn = seq.getAnnotation();
//                for (Iterator i = seqAn.keys().iterator(); i.hasNext(); ) {
//                    Object key = i.next();
//                    Object value = seqAn.getProperty(key);
//                    logger.info(key.toString() + ": " + value.toString());
//                }                
	}
 
	LinkedHashMap<String, ProteinSequence> protSequences = GenbankReaderHelper.readGenbankProteinSequence(protFile);
	for (ProteinSequence sequence : protSequences.values()) {
		logger.info( sequence.getSequenceAsString() );
	}
	/*
	 * Method 3: With the GenbankReader Object 
	 */		
	//Try reading with the GanbankReader
	FileInputStream is = new FileInputStream(dnaFile);
	GenbankReader<DNASequence, NucleotideCompound> dnaReader = new GenbankReader<DNASequence, NucleotideCompound>(
	        is, 
	        new GenericGenbankHeaderParser<DNASequence,NucleotideCompound>(),
	        new DNASequenceCreator(DNACompoundSet.getDNACompoundSet())
	);
	dnaSequences = dnaReader.process();
	is.close();
	logger.info(dnaSequences);
 
	is = new FileInputStream(protFile);
	GenbankReader<ProteinSequence, AminoAcidCompound> protReader = new GenbankReader<ProteinSequence, AminoAcidCompound>(
	        is,
	        new GenericGenbankHeaderParser<ProteinSequence,AminoAcidCompound>(),
	        new ProteinSequenceCreator(AminoAcidCompoundSet.getAminoAcidCompoundSet())
	);
	protSequences = protReader.process();
	is.close();
	logger.info(protSequences);
    } 
}
