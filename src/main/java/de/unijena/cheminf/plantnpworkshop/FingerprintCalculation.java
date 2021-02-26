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
import org.openscience.cdk.fingerprint.IBitFingerprint;
import org.openscience.cdk.fingerprint.PubchemFingerprinter;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
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
        //loading SD file from resources
        ClassLoader tmpClassLoader = DescriptorCalculation.class.getClassLoader();
        File tmpSDFile = new File(tmpClassLoader.getResource("COCONUTset-10.sdf").getFile());
        IChemObjectBuilder tmpBuilder = DefaultChemObjectBuilder.getInstance();
        IteratingSDFReader tmpReader = new IteratingSDFReader(new FileInputStream(tmpSDFile), tmpBuilder, true);

        //iterating molecules in file
        while(tmpReader.hasNext()) {
            IAtomContainer tmpMolecule = tmpReader.next();
            //reading attributes stored in SD, they are added as properties to the atom container automatically
            String tmpCOCONUTID = tmpMolecule.getProperty("COCONUT_ID");
            String tmpName = tmpMolecule.getProperty("Name");
            System.out.println("\n" + tmpName + " (" + tmpCOCONUTID + ")");

            //Preprocessing for PubChem fingerprint calculation (aromaticity detected, atom types configured, Hs explicit)
            AtomContainerManipulator.convertImplicitToExplicitHydrogens(tmpMolecule);
            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(tmpMolecule);
            Aromaticity tmpAromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.or(Cycles.all(), Cycles.cdkAromaticSet()));
            tmpAromaticity.apply(tmpMolecule);

            //PubChem fingerprint calculation (881 bits)
            PubchemFingerprinter tmpPubchemFingerprinter = new PubchemFingerprinter(tmpBuilder);
            IBitFingerprint tmpPubchemFingerprint = tmpPubchemFingerprinter.getBitFingerprint(tmpMolecule);
            //cardinality() returns the number of bits set to true in the fingerprint.
            System.out.println(tmpPubchemFingerprint.cardinality());
            //prints indices of the positive bits
            System.out.println(tmpPubchemFingerprint.asBitSet().toString());

            //TODO: next up, ECFP (preprocessing needed similar to PubChem FP?)
        }
    }

}
