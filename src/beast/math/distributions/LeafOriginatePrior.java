package beast.math.distributions;

import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.core.util.Log;
import beast.core.Input.Validate;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

@Description("Distribution on the originate of one or more leafs in a tree")
public class LeafOriginatePrior extends Distribution {
	public Input<Tree> treeInput = new Input<>("tree", "tree for which to use the leafs from", Validate.REQUIRED);
	public Input<String> valueInput = new Input<>("value", "tab delimited set of taxon-value pairs "
			+ "e.g. taxon1=100,taxon2=47,... where the taxa names are from the taxa in the tree and"
			+ "values indicate lower values of the age of the originate of the taxa", Validate.REQUIRED);

	Tree tree;
	Map<String, Double> map;
	
	final static double penalty = 1000;
	
	
	@Override
	public void initAndValidate() {		
		super.initAndValidate();
		tree = treeInput.get();
		
		Set<String> taxa = new LinkedHashSet<>();
		for(String taxon : tree.getTaxaNames()) {
			taxa.add(taxon);
		}
		
		map = new LinkedHashMap<>();
		String value = valueInput.get();
		String [] strs = value.split(",");
		for (String str : strs) {
			String [] strs2 = str.split("=");
			if (strs2.length > 1) {
				String taxon = strs2[0].trim();
				if (taxa.contains(taxon)) {
					String lowerBound = strs2[1].trim();
					Double d = Double.parseDouble(lowerBound);
					if (d > 0) {
						map.put(taxon, d);
					}
					taxa.remove(taxon);
				} else {
					Log.warning("LeafOriginatePrior has extra taxon in values: " + taxon);
				}
			}
		}
		
		for (String taxon : taxa) {
			Log.warning("LeafOriginatePrior has no value for: " + taxon);
		}
	}
	
	
	@Override
	public double calculateLogP() {
		logP = 0;
		for (Node leaf : tree.getExternalNodes()) {
			String taxon = leaf.getID();
			if (map.containsKey(taxon)) {
				double lowerBound = map.get(taxon);
				double h = leaf.getParent().getHeight();
				if (h < lowerBound) {
	                logP += -penalty * (lowerBound-h) * (lowerBound-h);
				}
			}
		}
		return logP;
	}
	
	
	@Override
	public List<String> getArguments() {return null;}

	@Override
	public List<String> getConditions() {return null;}

	@Override
	public void sample(State state, Random random) {}

}
