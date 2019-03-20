package beast.evolution.tree;

import java.util.*;

import org.apache.commons.math3.util.FastMath;

import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.State;
import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import beast.evolution.alignment.Alignment;
import beast.evolution.operators.MonoCladesMapping;
import beast.math.distributions.MRCAPrior;
import beast.math.distributions.MultiMRCAPriors;
import beast.math.distributions.MultiMonophyleticConstraint;
import test.beast.BEASTTestCase;


@Description("Tree prior that assumes a pure birth process with rate differences at specific clades (nested clades are not allowed)")
public class MultiYuleModel extends Distribution {
	public Input<Double> rhoInput = new Input<>("rho", "Extant sampling proportion, default 1", 1.0);

	public Input<RealParameter> lambdaInput = new Input<>("birthDiffRate", "birth rate, one for each clade", Validate.REQUIRED);
	public Input<RealParameter> gammaInput = new Input<>("gamma", "controls rate shifts (== the rate shift rate)", Validate.REQUIRED);

	public Input<MultiMonophyleticConstraint> constraintInput = new Input<>("constraint", "", Validate.REQUIRED);
	public Input<List<MRCAPrior>> excludeInput = new Input<>("exclude","list of clades which are included in constraints which are ignored - they " +
			"take the rate of the clade above.",
			new ArrayList<>());
	public Input<List<MRCAPrior>> includeInput = new Input<>("include","list of clades with their own rate. That is, all the branches from the " +
			"clade root to the top of other Monophletic clades in the constraints.",
			new ArrayList<>());

	final public Input<Boolean> excludeRootsInput = new Input<>("excludeRoots", "Roots of clades get the rate of their enclosing clade (default " +
			"false)", false);

	protected double rho;
	
	RealParameter gamma;
	RealParameter lambda;
	int[] nodeToCladeGroup = null;
	
	double [] oldLength;
	double [] oldRate;
	double [] oldLength1;
	double [] oldRate1;
	double [] logPContributionCache;
	double oldGamma = -1;

	Tree tree;
	MultiMonophyleticConstraint constraint;
	boolean excludeRoots;
	int[] roots = null;

	public void initAndValidate() {
		constraint = constraintInput.get();
		tree = constraint.treeInput.get();
		rho = rhoInput.get();
		gamma = gammaInput.get();
		lambda = lambdaInput.get();
		// for restore: init will resize to right length
		lambda.setDimension(tree.getNodeCount()-1);
		excludeRoots = excludeRootsInput.get();

//		// remove nesting clades from newick
//		String newick = newickInput.get();
//		while (newick.matches(".*\\(([^\\)]*)\\(.*")) {
//			newick = newick.replaceAll("\\(([^\\)]*)\\(", "($1");
//		}
//		while (newick.matches(".*\\)([^\\(]*)\\).*")) {
//			newick = newick.replaceAll("\\)([^\\(]*)\\)", "$1)");
//		}
//		newick = newick.replaceAll(";", "");
//		newick = "(" + newick + ");";
//		newickInput.setValue(newick, this);
		
		// checkValidity = checkValidityInput.get();
		// toThePresent = toThePresentInput.get();

		// rateShiftModel = rateShiftModelInput.get();
		//super.initAndValidate();
	}

	int getGroupIndex(Node node) {
		int nr = node.getNr();
		if( excludeRoots && ! node.isRoot() ) {
			if( roots[nodeToCladeGroup[nr]] == nr ) {
				nr = node.getParent().getNr();
			}
		}
		return nodeToCladeGroup[nr];
	}

	@Override
	public double calculateLogP() {
		
		if (nodeToCladeGroup == null) {
			initialise();
		}

		logP = 0;
		final double gamma = this.gamma.getValue();
		if( excludeRoots ) {
			for(int i = 0; i < roots.length; ++i) {
				if( roots[i] == -1 ) {
					continue;
				}
				Node n = tree.getNode(roots[i]);
				assert nodeToCladeGroup[roots[i]] == i;
				while( nodeToCladeGroup[n.getParent().getNr()] == i ) {
					n = n.getParent();
				}
				roots[i] = n.getNr();
				assert nodeToCladeGroup[roots[i]] == i;
			}
		}
		
		// gamma changed, need to recalculate the whole cache
		for (Node node : tree.getNodesAsArray()) {
			final int type = getGroupIndex(node); //nodeToCladeGroup[node.getNr()];
			final double lambda = this.lambda.getValue(type);
			final double t = node.getHeight();
			if (node.isRoot()) {
				if( rho != 1 ) {
					logP += -2 * FastMath.log1p(-p0(lambda, gamma, t));
				}
			} else {
				final double tp = node.getParent().getHeight();
				logP += logftip(lambda, gamma, tp, t);
				if( ! node.isLeaf() ) {
					// can cache some for speed later
					logP += FastMath.log(lambda);
				}
			}
		}

		if( false ) {
			logP = 0;
			for (Node node : tree.getNodesAsArray()) {
				final int type = getGroupIndex(node); //nodeToCladeGroup[node.getNr()];
				final double lambda = this.lambda.getValue(type);
				final double t = node.getHeight();
				final int i = node.getNr();
				if( node.isRoot() ) {
					logPContributionCache[i] = -2 * FastMath.log(1 - p0(lambda, gamma, t));
				} else {
					final double tp = node.getParent().getHeight();
					if( true ) {
						// test code
						logPContributionCache[i] = logftip(lambda, gamma, tp) - logftip(lambda, gamma, t);
						double tmp = logftip(lambda, gamma, tp, t);
						assert (Math.abs(logPContributionCache[i] - tmp)) < 1e-13;
					} else {
						logPContributionCache[i] = logftip(lambda, gamma, tp, t);
					}
				}
				oldLength[i] = t;
				oldRate[i] = lambda;
				logP += logPContributionCache[i];
			}
		}

		
/*		
		if (gamma != oldGamma) {
			// gamma changed, need to recalculate the whole cache
			for (Node node : tree.getNodesAsArray()) {
				final int type = getGroupIndex(node); //nodeToCladeGroup[node.getNr()];
				final double lambda = this.lambda.getValue(type);
				final double t = node.getHeight();
				final int i = node.getNr();
				if (node.isRoot()) {
//					logPContributionCache[i] = -2 * FastMath.log(1 - p0(lambda, gamma, t)) + 2 * logftip(lambda, gamma, t);
					logPContributionCache[i] = -2 * FastMath.log(1 - p0(lambda, gamma, t));
//				} else if (node.isLeaf()) {
//					logPContributionCache[i] = -logftip(lambda, gamma, t);
				} else {
//					logPContributionCache[i] = logftip(lambda, gamma, t);
					final double tp = node.getParent().getHeight();
					logPContributionCache[i] = logftip(lambda, gamma, tp) - logftip(lambda, gamma, t);
				}
				oldLength[i] = t;
				oldRate[i] = lambda;
				logP += logPContributionCache[i];
			}
			
		} else {
			// gamma unchanged, recalculate cache only if 
			// something changed
		
			for (Node node : tree.getNodesAsArray()) {
				final int type = getGroupIndex(node); //final int type = nodeToCladeGroup[node.getNr()];
				final double lambda = this.lambda.getValue(type);
				final double t = node.getHeight();
				final int i = node.getNr();
	
				if (oldLength[i] != t || oldRate[i] != lambda) {
					if (node.isRoot()) {
//						logPContributionCache[i] = -2 * FastMath.log(1 - p0(lambda, gamma, t)) + 2 * logftip(lambda, gamma, t);
						logPContributionCache[i] = -2 * FastMath.log(1 - p0(lambda, gamma, t));
//					} else if (node.isLeaf()) {
//						logPContributionCache[i] = -logftip(lambda, gamma, t);
					} else {
//						logPContributionCache[i] = logftip(lambda, gamma, t);
						final double tp = node.getParent().getHeight();
						logPContributionCache[i] = logftip(lambda, gamma, tp) - logftip(lambda, gamma, t);
					}
					oldLength[i] = t;
					oldRate[i] = lambda;
				}
				logP += logPContributionCache[i];
			}
		}
*/		
		for (Node node : tree.getNodesAsArray()) {
			if (node.getLength() < 1e-6 && !node.isRoot() && node.getLength() > 0) {
				logP -= 1.0/node.getLength(); 
			}
		}
		oldGamma = gamma;
		
//		double logP2 = calculateLogP2();
//		if (Math.abs(logP - logP2) > 1e-9) {
//			System.err.print("x");
//		}
		
		return logP;
	}
	
//	public double calculateLogP2() throws Exception {
//
//		if (nodeToCladeGroup == null) {
//			initialise();
//		}
//
//		double logP = 0;
//		double gamma = this.gamma.getValue();
//
//		for (Node node : tree.getNodesAsArray()) {
//			int type = nodeToCladeGroup[node.getNr()];
//			double lambda = this.lambda.getValue(type);
//			double t = node.getHeight();
//			int i = node.getNr();
//			if (node.isRoot()) {
//				logP += - 2 * FastMath.log(1 - p0(lambda, gamma, t));
//			} else {
//				Node parent = node.getParent();
//				double t1 = parent.getHeight();
//				int type1 = nodeToCladeGroup[parent.getNr()];
//				double lambda1 = this.lambda.getValue(type1);
//				logP += logftip(lambda1, gamma, t1) - logftip(lambda, gamma, t);
//			}
//		}
//		if (true) {
//			return logP;
//		}
//
//		if (gamma != oldGamma) {
//			// gamma changed, need to recalculate the whole cache
//			for (Node node : tree.getNodesAsArray()) {
//				int type = nodeToCladeGroup[node.getNr()];
//				double lambda = this.lambda.getValue(type);
//				double t = node.getHeight();
//				int i = node.getNr();
//				if (node.isRoot()) {
//					logPContributionCache[i] = - 2 * FastMath.log(1 - p0(lambda, gamma, t));
//				} else {
//					Node parent = node.getParent();
//					double t1 = parent.getHeight();
//					int type1 = nodeToCladeGroup[parent.getNr()];
//					double lambda1 = this.lambda.getValue(type1);
//					logPContributionCache[i] = logftip(lambda1, gamma, t1) - logftip(lambda, gamma, t);
//					oldLength1[i] = t1;
//					oldRate1[i] = lambda1;
//				}
//				oldLength[i] = t;
//				oldRate[i] = lambda;
//				logP += logPContributionCache[i];
//			}
//
//		} else {
//			// gamma unchanged, recalculate cache only if
//			// something changed
//
//			for (Node node : tree.getNodesAsArray()) {
//				int type = nodeToCladeGroup[node.getNr()];
//				double lambda = this.lambda.getValue(type);
//				double t = node.getHeight();
//				int i = node.getNr();
//
//				if (node.isRoot()) {
//					if (oldLength[i] != t || oldRate[i] != lambda) {
//						logPContributionCache[i] = - 2 * FastMath.log(1 - p0(lambda, gamma, t));
//						oldLength[i] = t;
//						oldRate[i] = lambda;
//					}
//				} else {
//					Node parent = node.getParent();
//					double t1 = parent.getHeight();
//					int type1 = nodeToCladeGroup[parent.getNr()];
//					double lambda1 = this.lambda.getValue(type1);
//					if (oldLength[i] != t || oldRate[i] != lambda || oldLength1[i] != t1 || oldRate1[i] != lambda1) {
//						logPContributionCache[i] = logftip(lambda1, gamma, t1) - logftip(lambda, gamma, t);
//						oldLength1[i] = t1;
//						oldRate1[i] = lambda1;
//						oldLength[i] = t;
//						oldRate[i] = lambda;
//					}
//				}
//				logP += logPContributionCache[i];
//			}
//		}
//		oldGamma = gamma;
//
//		return logP;
//	}

	private void adjust(Node n, int masterGroup, boolean split, Set<Integer> ownRate, Set<Integer> ignored) {
		boolean inherit = true;
		int nid = n.getNr();
		if( ! n.isLeaf() ) {
			// Non-leaf case
			final int ng = nodeToCladeGroup[nid];
			final int png = n.isRoot() ? ng : nodeToCladeGroup[n.getParent().getNr()];
			if( ng != png ) {
				if( ignored.contains(nid) ) {
					// inherit = true;
				} else if( ownRate.contains(nid) ) {
					inherit = false;
					masterGroup = ng;
					split = true;
				} else {
					if( split ) {
						inherit = false;
						// own rate by default
						masterGroup = ng;
						split = false;
					}
				}
			} else {
				assert ! ignored.contains(nid) &&! ownRate.contains(nid);
			}

			for( Node c : n.getChildren() ) {
				adjust(c, masterGroup, split, ownRate, ignored);
			}
		} else {
			// Leaf case
			final int ng = nodeToCladeGroup[nid];
			if (ng != -1 && constraint.getConstraints().get(ng).size() == 1) {
	                    // This is a leaf which belongs to a monophyly constraint
			    // of size 1.  This is almost certainly an isolate, which
			    // we'd like to get its own rate.
			    inherit = false;
		        }
		}
		if( inherit ) {
			nodeToCladeGroup[nid] = masterGroup;
		}
	}

	private void adjustCladeMapping() {
		Set<Integer> ownRate = new HashSet<>();
		Set<Integer> ignored = new HashSet<>();

		for (MRCAPrior prior : excludeInput.get()) {
			int nodeNr = prior.getCommonAncestor().getNr();
			ignored.add(nodeNr);
		}
		for (MRCAPrior prior : includeInput.get()) {
			int nodeNr = prior.getCommonAncestor().getNr();
			ownRate.add(nodeNr);
		}

		final Node root = tree.getRoot();
		adjust(root, nodeToCladeGroup[root.getNr()], true, ownRate, ignored);
	}

	void initialise() {
		if (Double.isInfinite(constraint.calculateLogP())) {
			throw new RuntimeException("Starting tree does not conform to monopnyletic constraints");
		}
		nodeToCladeGroup = MonoCladesMapping.setupNodeGroup(tree, constraint);
		adjustCladeMapping();

		// Renumber the clade groups from 0 to n-1
		// First, make sure that the root birthrate is numbered 0
		Map<Integer,Integer> map = new HashMap<>();
		int dim = 0;
		final int rootGroup = nodeToCladeGroup[tree.getRoot().getNr()];
		map.put(rootGroup, dim);
		dim++;
		// Then renumber the other groups in a random order
		for (int iNode2Group : nodeToCladeGroup) {
			if( !map.containsKey(iNode2Group) ) {
				map.put(iNode2Group, dim);
				dim++;
			}
		}
		for (int i = 0; i < nodeToCladeGroup.length; i++) {
			nodeToCladeGroup[i] = map.get(nodeToCladeGroup[i]);
		}			

		lambda.setDimension(dim);

		if( excludeRoots ) {
			roots = new int[dim];
			for (int i = 0; i < nodeToCladeGroup.length; i++) {
				Node n = tree.getNode(i);
				if( n.isRoot() ) {
					roots[ nodeToCladeGroup[i] ] = -1;
				} else {
					if( nodeToCladeGroup[i] != nodeToCladeGroup[n.getParent().getNr()] ) {
						if (!n.isLeaf()) {
							roots[nodeToCladeGroup[i]] = i;
						} else {
							// leaf nodes should keep their own rate
							roots[nodeToCladeGroup[i]] = -1;
						}
					}
				}
			}
		}
		oldLength = new double[tree.getNodeCount()];
		oldRate = new double[tree.getNodeCount()];
		oldLength1 = new double[tree.getNodeCount()];
		oldRate1 = new double[tree.getNodeCount()];
		logPContributionCache = new double[tree.getNodeCount()];
		
		if (gamma.getLower() < 0) {
			gamma.setLower(0.0);
		}
		if (gamma.getUpper() > 1.0) {
			gamma.setUpper(1.0);
		}
		
		if (constraint instanceof MultiMRCAPriors) {
			String [] birthRateID = new String[dim];
			for (MRCAPrior prior : ((MultiMRCAPriors) constraint).calibrationsInput.get()) {
				if (excludeInput.get().indexOf(prior) == -1 ) {
					int k = prior.getCommonAncestor().getNr();
					if (nodeToCladeGroup[k] != -1) {
						birthRateID[nodeToCladeGroup[k]] = prior.getID().replaceAll(".prior", "");
					}
				}
			}
			for (int i = 0; i < dim; i++) {
				int tot = 0;
				for (int aNodeToCladeGroup : nodeToCladeGroup) {
					tot += (aNodeToCladeGroup == i) ? 1 : 0;
				}
				Log.warning("birth rate " + (i+1) + " => " + (birthRateID[i] == null ? "[root clade]" : birthRateID[i]) + "(" + tot + " nodes)");
			}
		}
	}

	protected double p0(double lambda, /* double mu, */double gamma, double time) {
		//double c = Math.sqrt(Math.pow(gamma + 0 + lambda, 2) - 4 * 0 * lambda * (1 - 0));
		//double c = gamma + lambda;
		//double x1 = -(gamma + 0 + lambda + c) / 2;
		double x1 = -(gamma + lambda);
		//double x2 = -(gamma + 0 + lambda - c) / 2; == 0
		//double A1 = lambda * (1 - rho) + x1;
		double A1 = -lambda * rho - gamma;
		double A2 = lambda * (1 - rho);
		//double res = (A1 * 0 - A2 * x1 * Math.exp(-c * time)) / (lambda * (A2 * Math.exp(-c * time) - A1));
		double res = (A2 * x1 * FastMath.exp(x1 * time)) / (lambda * (A2 * FastMath.exp(x1 * time) - A1));
		return res;
	}

	protected double logftip(double lambda, /* double mu, */double gamma, double time) {
		//double c = Math.sqrt(Math.pow(gamma + 0 + lambda, 2) - 4 * 0 * lambda * (1 - 0));
		//double c = gamma + lambda;
		double x1 = -(gamma + lambda);
		//double x2 = -(gamma + 0 + lambda - c) / 2; == 0
		//double A1 = lambda * (1 - rho) + x1;
		double A1 = -lambda * rho - gamma;
		double A2 = lambda * (1 - rho);
		double y = A2 * FastMath.exp(x1 * time) - A1;
		double res = x1 * x1 / (FastMath.exp(-x1 * time) * y * y);
		return FastMath.log(res);
	}

	protected double logftip(double lambda, /* double mu, */double gamma, double time1, double time2) {
		double v;
		final double x1 = -(gamma + lambda);
		if( rho == 1 ) {
			v = x1 * (time1 - time2);
		} else {
			final double mlamrho = -lambda * rho;
			final double A1 = mlamrho - gamma;
			final double A2 = lambda + mlamrho;

			final double x1t1 = x1 * time1;
			final double y1 = A2 * FastMath.exp(x1t1) - A1;

			final double x1t2 = x1 * time2;
			final double y2 = A2 * FastMath.exp(x1t2) - A1;

			v = (x1t1 - x1t2) + 2 * FastMath.log(y2 / y1);
		}
		return v;
	}
//	private double logPNode(Node node, double tx) {
//		double lambda, /* mu, */ti, te;
//		int type;
//		double gamma = this.gamma.getValue();
//		double logP = 0;
//
//		type = nodeToCladeGroup[node.getNr()];
//		lambda = this.lambda.getValue(type);
//		te = node.getHeight();
//		ti = node.getParent().getHeight();
//		logP += logftip(lambda, gamma, ti) - logftip(lambda, gamma, te);
//
//		if (!node.isLeaf())
//			logP += logPNode(node.getLeft(), tx) + logPNode(node.getRight(), tx);
//		return (logP);
//	}

	@Override
	public List<String> getArguments() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public List<String> getConditions() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void sample(State state, Random random) {
		// TODO Auto-generated method stub

	}
	
	
	@Override
	protected boolean requiresRecalculation() {
		return true;
	}

	
	public static void main(String[] args) throws Exception {
		Alignment data = BEASTTestCase.getAlignment();
	    Tree tree = BEASTTestCase.getTree(data);
	    
		MultiYuleModel myd = new MultiYuleModel();
		myd.initByName("tree", tree, 
				"newick", "(human:0.024003,chimp:0.010772,bonobo:0.010772),gorilla:0.036038,orangutan:0.069125,siamang:0.099582;",
				"birthDiffRate", "0.1",
				"gamma", "0.5");
		
		System.err.println("logP = " + myd.calculateLogP());
	}
}
