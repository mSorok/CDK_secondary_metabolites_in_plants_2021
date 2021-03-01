/*
 * MIT License
 *
 * Copyright (c) 2021. Maria Sorokina, Aziz M. Yirik, Jonas Schaub, Christoph Steinbeck
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

package de.unijena.cheminf.plantnpworkshop;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fingerprint.CircularFingerprinter;
import org.openscience.cdk.fingerprint.IBitFingerprint;
import org.openscience.cdk.fingerprint.PubchemFingerprinter;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.similarity.Tanimoto;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;

/**
 * Demo class for calculating PubChem and EC fingerprints using CDK.
 *
 * @author Jonas Schaub
 */
public class FingerprintCalculation {

    /**
     * COCONUT subset molecules are loaded from SD file and their PubChem and Extended Connectivity fingerprints
     * calculated and reported on console.
     *
     * @param args the command line arguments (none required)
     */
    public static void main(String[] args) throws FileNotFoundException, CDKException {

        //*loading SD file from resources*
        File tmpSDFile = new File("src/main/resources/COCONUTset-10.sdf");
        //a factory class to provide implementation independent ICDKObjects, needed for SDF parsing
        IChemObjectBuilder tmpBuilder = DefaultChemObjectBuilder.getInstance();
        //skip: true -> erroneous entries will be skipped
        IteratingSDFReader tmpReader = new IteratingSDFReader(new FileInputStream(tmpSDFile), tmpBuilder, true);

        //*iterating molecules in file*
        while(tmpReader.hasNext()) {
            IAtomContainer tmpMolecule = tmpReader.next();
            //reading attributes of molecules stored in SDF, they are added as properties to the atom container automatically
            String tmpCOCONUTID = tmpMolecule.getProperty("COCONUT_ID");
            String tmpName = tmpMolecule.getProperty("Name");
            System.out.println("\n\n" + tmpName + " (" + tmpCOCONUTID + ")");

            //*PubChem fingerprint calculation*
            //preprocessing required: Hs explicit, atom types configured, aromaticity detected
            SmilesGenerator tmpSmiGen = new SmilesGenerator(SmiFlavor.Unique | SmiFlavor.UseAromaticSymbols);
            System.out.println("\tSMILES representation before preprocessing: " + tmpSmiGen.create(tmpMolecule));
            //aromaticity model is constructed from electron donation model and cycle finder
            AtomContainerManipulator.convertImplicitToExplicitHydrogens(tmpMolecule);
            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(tmpMolecule);
            Aromaticity tmpAromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.cdkAromaticSet());
            tmpAromaticity.apply(tmpMolecule);
            System.out.println("\tSMILES representation AFTER preprocessing: " + tmpSmiGen.create(tmpMolecule));
            PubchemFingerprinter tmpPubChemFingerprinter = new PubchemFingerprinter(tmpBuilder);
            IBitFingerprint tmpPubChemFingerprint = tmpPubChemFingerprinter.getBitFingerprint(tmpMolecule);
            System.out.println("\n\tPubChem fingerprint:");
            System.out.println("\t\tNumber of positive bits: " + tmpPubChemFingerprint.cardinality());
            System.out.println("\t\tIndices of positive bits: " + tmpPubChemFingerprint.asBitSet().toString());

            //*ECFP calculation*
            //Circular fingerprints: for generating fingerprints that are functionally equivalent to ECFP-2/4/6 and FCFP-2/4/6 fingerprints
            //Default constructor: uses the ECFP6 type.
            CircularFingerprinter tmpECFPrinter = new CircularFingerprinter();
            //aromaticity detection and atom typing is done internally
            //implicit vs. explicit hydrogens are handled, i.e. it doesn't matter whether the incoming molecule is hydrogen suppressed or not.
            //Calculates the circular fingerprint for the given IAtomContainer, and folds the result into a single bitset (see getSize()).
            IBitFingerprint tmpECFPrint = tmpECFPrinter.getBitFingerprint(tmpMolecule);
            System.out.println("\n\tCircular fingerprint (ECFP6):");
            System.out.println("\t\tNumber of positive bits: " + tmpECFPrint.cardinality());
            System.out.println("\t\tIndices of positive bits: " + tmpECFPrint.asBitSet().toString());

            //*Tanimoto calculation*
            SmilesParser tmpSmiPar = new SmilesParser(tmpBuilder);
            //SMILES representation of Flower Of Paradise (CNP0218319)
            IAtomContainer tmpFlowerOfParadise = tmpSmiPar.parseSmiles("O=C1C=C(O)C(=O)C=2C=CC=CC12");
            IBitFingerprint tmpFOPECFPrint = tmpECFPrinter.getBitFingerprint(tmpFlowerOfParadise);
            AtomContainerManipulator.convertImplicitToExplicitHydrogens(tmpFlowerOfParadise);
            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(tmpFlowerOfParadise);
            tmpAromaticity.apply(tmpFlowerOfParadise);
            IBitFingerprint tmpFOPPubChemFingerprint = tmpPubChemFingerprinter.getBitFingerprint(tmpFlowerOfParadise);
            double tmpPubChemTanimoto = Tanimoto.calculate(tmpPubChemFingerprint, tmpFOPPubChemFingerprint);
            System.out.println("\n\tTanimoto similarity (using PubChem FP) to the Flower of Paradise: " + String.format("%,.2f", tmpPubChemTanimoto));
            double tmpECFPTanimoto = Tanimoto.calculate(tmpECFPrint, tmpFOPECFPrint);
            System.out.println("\tTanimoto similarity (using ECFP) to the Flower of Paradise: " + String.format("%,.2f", tmpECFPTanimoto));
        } //end of while()
    } //end of main()
}
