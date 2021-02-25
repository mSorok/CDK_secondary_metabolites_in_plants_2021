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
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.qsar.DescriptorValue;
import org.openscience.cdk.qsar.descriptors.molecular.ALOGPDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.PetitjeanNumberDescriptor;
import org.openscience.cdk.qsar.result.DoubleArrayResult;
import org.openscience.cdk.qsar.result.DoubleResult;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.util.Objects;

/**
 *
 */
public class DescriptorCalculation {
    /**
     *
     * @param args
     */
    public static void main(String[] args) throws FileNotFoundException, CDKException {
        ClassLoader tmpClassLoader = DescriptorCalculation.class.getClassLoader();
        File tmpSDFile = new File(tmpClassLoader.getResource("COCONUTset-10.sdf").getFile());
        IChemObjectBuilder tmpBuilder = DefaultChemObjectBuilder.getInstance();
        IteratingSDFReader tmpReader = new IteratingSDFReader(new FileInputStream(tmpSDFile), tmpBuilder, true);
        Aromaticity tmpAromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.or(Cycles.all(), Cycles.cdkAromaticSet()));
        ALOGPDescriptor tmpALogPDescriptor = new ALOGPDescriptor();
        tmpALogPDescriptor.initialise(tmpBuilder);
        PetitjeanNumberDescriptor tmpPetitjeanNumberDescriptor = new PetitjeanNumberDescriptor();
        tmpPetitjeanNumberDescriptor.initialise(tmpBuilder);
        String tmpCOCONUTID;
        String tmpName;
        IAtomContainer tmpMolecule;
        while(tmpReader.hasNext()) {
            tmpMolecule = tmpReader.next();
            if (Objects.isNull(tmpMolecule)) {
                break;
            }
            tmpCOCONUTID = tmpMolecule.getProperty("COCONUT_ID");
            tmpName = tmpMolecule.getProperty("Name");
            System.out.println("\n" + tmpName + " (" + tmpCOCONUTID + ")");
            //for alogp calculation, Hs must be explicit and aromaticity detected
            AtomContainerManipulator.convertImplicitToExplicitHydrogens(tmpMolecule);
            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(tmpMolecule);
            tmpAromaticity.apply(tmpMolecule);
            //System.out.println((new SmilesGenerator(SmiFlavor.Unique | SmiFlavor.UseAromaticSymbols)).create(tmpMolecule));
            //calculates 3 values, ALogP - Ghose-Crippen LogKow, ALogP2, amr - molar refractivity
            //result type is an array of doubles
            DescriptorValue tmpALogPDescriptorResults = tmpALogPDescriptor.calculate(tmpMolecule);
            String[] tmpALogPDescriptorNames = tmpALogPDescriptorResults.getNames();
            DoubleArrayResult tmpALogPDescriptorResultsValues = (DoubleArrayResult) tmpALogPDescriptorResults.getValue();
            for (int i = 0; i < tmpALogPDescriptorNames.length; i++) {
                System.out.println("\t" + tmpALogPDescriptorNames[i] + ": " + String.format("%,.2f", tmpALogPDescriptorResultsValues.get(i)));
            }
            DescriptorValue tmpPetitjeanNumberDescriptorResult = tmpPetitjeanNumberDescriptor.calculate(tmpMolecule);
            DoubleResult tmpPetitjeanNumberDescriptorResultValue = (DoubleResult) tmpPetitjeanNumberDescriptorResult.getValue();
            System.out.println("\t" + tmpPetitjeanNumberDescriptorResult.getNames()[0] + ": " + String.format("%,.2f", tmpPetitjeanNumberDescriptorResultValue.doubleValue()));
            //https://github.com/mSorok/COCONUT/blob/master/src/main/java/de/unijena/cheminf/npopensourcecollector/services/MolecularFeaturesComputationService.java
        }
    }
}
