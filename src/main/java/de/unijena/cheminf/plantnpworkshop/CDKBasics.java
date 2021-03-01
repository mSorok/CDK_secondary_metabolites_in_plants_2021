package de.unijena.cheminf.plantnpworkshop;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;

import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.depict.DepictionGenerator;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.MDLV3000Reader;
import org.openscience.cdk.io.SDFWriter;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;

public class CDKBasics {

	public static SDFWriter outputFile;
	public static int count;
	
	/**
	 * SMILES parser and generator functions.
	 */
	
	/**
	 * Building an IAtomContainer from a given SMILES.
	 * 
	 * @param smiles String SMILES representation of a molecule
	 * @return IAtomContainer 
	 * @throws InvalidSmilesException
	 */
	
	public static IAtomContainer smiles2IAtomContainer(String smiles) throws InvalidSmilesException {
		SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
		return sp.parseSmiles(smiles);
	}
	
	/**
	 * For an IAtomContaineri getting the SMILES representation.
	 * 
	 * @param ac IAtomContainer molecule 
	 * @return String SMILES
	 * @throws CDKException
	 */
	
	public static String getSMILES(IAtomContainer ac) throws CDKException {
		
		/**
		 * There are different SMILES types. Herei we prefered 'Unique' SMILES.
		 */
		 
		SmilesGenerator sg = new SmilesGenerator(SmiFlavor.Unique);
		return sg.create(ac); 
	}
	
	/**
	 * Writing - Reading SDF Files
	 */
	
	/**
	 * Setting the directory path for the output SDF file.
	 * 
	 * @param path String directory path
	 * @throws IOException
	 */
	
	public static void setSDFfile(String path) throws IOException {
		outputFile = new SDFWriter(new FileWriter(path+"output.sdf"));
	}
	
	/**
	 * Reading a chemistry file, such as MOL and SDF, and building a List<IAtomContainer> for all the molecules.
	 * @param path String SDF file directory
	 * @return List<IAtomContainer>
	 * @throws IOException 
	 */
	
	public static List<IAtomContainer> chemFile2IAtomContainers(String path) throws IOException{
		List<IAtomContainer> list= new ArrayList<IAtomContainer>();
		File sdfFile = new File(path);
		IteratingSDFReader reader = new IteratingSDFReader(new FileInputStream(sdfFile), DefaultChemObjectBuilder.getInstance());
		while (reader.hasNext()) {
		   list.add((IAtomContainer)reader.next());
		}
		reader.close();
		return list;
	}
	
	/**
	 * Reading a MOL file with single molecule. 
	 * 
	 * @param path String directory path
	 * @return IAtomContainer
	 * @throws CDKException
	 * @throws IOException
	 */
	
	public static IAtomContainer readMOLfile(String path) throws CDKException, IOException {
		InputStream in = new FileInputStream(path);
        MDLV3000Reader reader = new MDLV3000Reader(in);
        IAtomContainer mol = (IAtomContainer) reader.read(new AtomContainer());
        reader.close();
        return mol;
	}
	
	/**
	 * Molecule depiction.
	 * 
	 * @param molecule IAtomContainer molecule
	 * @param path String directory for the output png file.
	 * @throws CloneNotSupportedException
	 * @throws CDKException
	 * @throws IOException
	 */
	
	public static void depict(IAtomContainer molecule, String path) throws CloneNotSupportedException, CDKException, IOException{
		DepictionGenerator depiction = new DepictionGenerator();
		
		/**
		 * Users can set different parameters and functions.
		 */
		depiction.withAtomColors().withCarbonSymbols().withSize(1000, 1000).withZoom(20).depict(molecule).writeTo(path+count+".png"); //Counting the number of read IAtomContainers.
	}
	
	public static void main(String[] arguments) throws IOException, CDKException {	
	
	}
}
