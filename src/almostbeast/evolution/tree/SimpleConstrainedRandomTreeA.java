package almostbeast.evolution.tree;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

import org.apache.commons.math.MathException;

import almostbeast.math.distributions.AlmostDistribution;
import beast.base.core.BEASTInterface;
import beast.base.core.Description;
import beast.base.core.Log;
import beast.base.evolution.alignment.TaxonSet;
import beastlabs.evolution.tree.SimpleConstrainedRandomTree;
import beast.base.evolution.tree.MRCAPrior;
import beast.base.evolution.tree.Node;
import beast.base.inference.distribution.ParametricDistribution;
import beast.base.util.Randomizer;

@Description("Like SimpleConstrainedRandomTree, but attempting to take Almost distributions in account")
public class SimpleConstrainedRandomTreeA extends SimpleConstrainedRandomTree {
	
	
	@Override
    public void doTheWork() {
        // find taxon sets we are dealing with
        taxonSets = new ArrayList<>();
        m_bounds = new ArrayList<>();
        distributions = new ArrayList<>();
        taxonSetIDs = new ArrayList<>();
        List<Boolean> onParent = new ArrayList<>();
        lastMonophyletic = 0;

        if (taxaInput.get() != null) {
            sTaxa.addAll(taxaInput.get().getTaxaNames());
        } else {
            sTaxa.addAll(m_taxonset.get().asStringList());
        }

        // pick up constraints from outputs, m_inititial input tree and output tree, if any
        List<MRCAPrior> calibrations = new ArrayList<MRCAPrior>();
        calibrations.addAll(calibrationsInput.get());

        // pick up constraints in m_initial tree
        for (final Object plugin : getOutputs()) {
            if (plugin instanceof MRCAPrior && !calibrations.contains(plugin) ) {
                calibrations.add((MRCAPrior) plugin);
            }
        }

        if (m_initial.get() != null) {
            for (final Object plugin : m_initial.get().getOutputs()) {
                if (plugin instanceof MRCAPrior && !calibrations.contains(plugin)) {
                    calibrations.add((MRCAPrior) plugin);
                }
            }
        }

        for (final MRCAPrior prior : calibrations) {
        	if (prior.getID() != null && prior.getID().startsWith("Afro-Asiatic.o")) {
        		int h = 3;
        		h++;
        	}
            final TaxonSet taxonSet = prior.taxonsetInput.get();
            if (taxonSet != null && !prior.onlyUseTipsInput.get()) {
	            final Set<String> bTaxa = new LinkedHashSet<>();
	        	if (taxonSet.asStringList() == null) {
	        		taxonSet.initAndValidate();
	        	}
	            for (final String sTaxonID : taxonSet.asStringList()) {

	                if (!sTaxa.contains(sTaxonID)) {
	                    throw new IllegalArgumentException("Taxon <" + sTaxonID + "> could not be found in list of taxa. Choose one of " +
                                Arrays.toString(sTaxa.toArray(new String[sTaxa.size()])));
	                }
	                bTaxa.add(sTaxonID);
	            }
	            final ParametricDistribution distr = prior.distInput.get();
	            final Bound bounds = new Bound();
	            if (distr != null) {
	        		List<BEASTInterface> plugins = new ArrayList<>();
	        		distr.getPredecessors(plugins);
	        		for (int i = plugins.size() - 1; i >= 0 ; i--) {
	        			plugins.get(i).initAndValidate();
	        		}
	                try {
	                	if (distr instanceof AlmostDistribution) {
	                        bounds.lower = ((AlmostDistribution)distr).getLowerTarget();
			                bounds.upper = ((AlmostDistribution)distr).getUpperTarget();
	                	} else {
	                        bounds.lower = Math.max(distr.inverseCumulativeProbability(0.0), 0.0);
			                bounds.upper = distr.inverseCumulativeProbability(1.0);
	                	}
                        assert bounds.lower <= bounds.upper;
					} catch (MathException e) {
						Log.warning.println("Could not set bounds in SimpleRandomTree::doTheWork : " + e.getMessage());
					}
	            }

	            if (prior.isMonophyleticInput.get() || bTaxa.size() == 1 ) {
	                // add any monophyletic constraint
                    boolean isDuplicate = false;
                    for(int k = 0; k < lastMonophyletic; ++k) {
                        // assert prior.useOriginateInput.get().equals(onParent.get(k)) == (prior.useOriginateInput.get() == onParent.get(k));
                        if( bTaxa.size() == taxonSets.get(k).size() && bTaxa.equals(taxonSets.get(k)) &&
                                prior.useOriginateInput.get().equals(onParent.get(k)) ) {
                            if( distr != null ) {
                                if( distributions.get(k) == null ) {
                                    distributions.set(k, distr);
                                    m_bounds.set(k, bounds);
                                    taxonSetIDs.set(k, prior.getID());
                                    System.out.println(prior.getID() + " " + bounds);
                                }
                            } else {
                                System.out.println("Duplicate " + prior.getID() + " (" + Arrays.toString(bTaxa.toArray()) + ") originate=" + prior.useOriginateInput.get());                            	
                            }
                            isDuplicate = true;
                        }
                    }
                    if( ! isDuplicate ) {
                        taxonSets.add(lastMonophyletic, bTaxa);
                        distributions.add(lastMonophyletic, distr);
                        onParent.add(lastMonophyletic, prior.useOriginateInput.get());
                        m_bounds.add(lastMonophyletic, bounds);
                        taxonSetIDs.add(lastMonophyletic, prior.getID());
                        lastMonophyletic++;
                    }
	            } else {
	                // only calibrations with finite bounds are added
	                if (!Double.isInfinite(bounds.lower) || !Double.isInfinite(bounds.upper)) {
	                    taxonSets.add(bTaxa);
	                    distributions.add(distr);
	                    m_bounds.add(bounds);
	                    taxonSetIDs.add(prior.getID());
                        onParent.add(prior.useOriginateInput.get());
	                }
	            }
            }
        }

        if( ICC ) {
            for (int i = 0; i < lastMonophyletic; i++) {
                final Set<String> ti = taxonSets.get(i);
                for (int j = i + 1; j < lastMonophyletic; j++) {
                    final Set<String> tj = taxonSets.get(j);
                    boolean i_in_j = tj.containsAll(ti);
                    boolean j_in_i = ti.containsAll(tj);
                    if( i_in_j || j_in_i ) {
                        boolean ok = true;
                        if( i_in_j && j_in_i ) {
                            ok = (boolean) (onParent.get(i)) != (boolean) onParent.get(j);
                        }
                        assert ok : "" + i + ' ' + j + ' ' + ' ' + taxonSetIDs.get(i) + ' ' + taxonSetIDs.get(j);
                    } else {
                        Set<String> tmp = new HashSet<>(tj);
                        tmp.retainAll(ti);
                        assert tmp.isEmpty();
                    }
                }
            }
        }

        // assume all calibration constraints are Monophyletic
        // TODO: verify that this is a reasonable assumption
        lastMonophyletic = taxonSets.size();

        // sort constraints in increasing set inclusion order, i.e. such that if taxon set i is subset of taxon set j, then i < j
//        for (int i = 0; i < lastMonophyletic; i++) {
//            for (int j = i + 1; j < lastMonophyletic; j++) {
//                final Set<String> taxai = taxonSets.get(i);
//                final Set<String> taxaj = taxonSets.get(j);
//                if (taxai.size() > taxaj.size() ||
//                	(taxai.size() == taxaj.size() && onParent.get(i))) {
//                    swap(taxonSets, i, j);
//                    swap(distributions, i, j);
//                    swap(m_bounds, i, j);
//                    swap(taxonSetIDs, i, j);
//                    swap(onParent, i, j);                	
//                }
//            }
//        }
        for (int i = 0; i < lastMonophyletic; i++) {
            for (int j = i + 1; j < lastMonophyletic; j++) {

                final Set<String> taxai = taxonSets.get(i);
                final Set<String> taxaj = taxonSets.get(j);
                Set<String> intersection = new LinkedHashSet<>(taxai);
                intersection.retainAll(taxaj);

                if (intersection.size() > 0) {
                    final boolean bIsSubset  = taxai.containsAll(taxaj);
                    final boolean bIsSubset2 = taxaj.containsAll(taxai);
                    // sanity check: make sure either
                    // o taxonset1 is subset of taxonset2 OR
                    // o taxonset1 is superset of taxonset2 OR
                    // o taxonset1 does not intersect taxonset2
                    if (!(bIsSubset || bIsSubset2)) {
                        throw new IllegalArgumentException("333: Don't know how to generate a Random Tree for taxon sets that intersect, " +
                                "but are not inclusive. Taxonset " + (taxonSetIDs.get(i) == null ? taxai :  taxonSetIDs.get(i)) + 
                                		" and " + (taxonSetIDs.get(j) == null ? taxaj : taxonSetIDs.get(j)));
                    }
                    // swap i & j if b1 subset of b2. If equal sub-sort on 'useOriginate'
                    if (bIsSubset && (!bIsSubset2 || (onParent.get(i) && !onParent.get(j)) ) ) {
                        swap(taxonSets, i, j);
                        swap(distributions, i, j);
                        swap(m_bounds, i, j);
                        swap(taxonSetIDs, i, j);
                        swap(onParent, i, j);
                    }
                }
            }
        }

        if( ICC ) {
            for (int i = 0; i < lastMonophyletic; i++) {
                final Set<String> ti = taxonSets.get(i);
                for (int j = i + 1; j < lastMonophyletic; j++) {
                    final Set<String> tj = taxonSets.get(j);
                    boolean ok = tj.containsAll(ti);
                    if( ok ) {
                        ok = !tj.equals(ti) || (!onParent.get(i) && onParent.get(j));
                        assert ok : "" + i + ' ' + j + ' ' + tj.equals(ti) + ' ' + taxonSetIDs.get(i) + ' ' + taxonSetIDs.get(j);
                    } else {
                        Set<String> tmp = new HashSet<>(tj);
                        tmp.retainAll(ti);
                        assert tmp.isEmpty();
                    }
                }
            }
        }

        for (int i = 0; i < lastMonophyletic; i++) {
            if( onParent.get(i) ) {
                // make sure it is after constraint on node itself, if such exists
                assert( ! (i + 1 < lastMonophyletic && taxonSets.get(i).equals(taxonSets.get(i + 1)) && onParent.get(i) && !onParent.get(i+1) ) );
                // find something to attach to ....
                // find enclosing clade, if any. pick a non-intersecting clade in the enclosed without an onParent constraint, or one whose
                // onParent constraint is overlapping.
                final Set<String> iTaxa = taxonSets.get(i);
                int j = i+1;
                Set<String> enclosingTaxa = sTaxa;
                {
                    String someTaxon = iTaxa.iterator().next();
                    for (/**/; j < lastMonophyletic; j++) {
                        if( taxonSets.get(j).contains(someTaxon) ) {
                            enclosingTaxa = taxonSets.get(j);
                            break;
                        }
                    }
                }
                if (i == 54) {
                	int h = 3;
                	h++;
                }
                final int enclosingIndex = (j == lastMonophyletic) ? j : j;
                Set<String> candidates = new HashSet<>(enclosingTaxa);
                candidates.removeAll(iTaxa);
                Set<Integer> candidateClades = new HashSet<>(5);
                List<String> canTaxa = new ArrayList<>();
                for( String c : candidates ) {
                    for(int k = enclosingIndex-1; k >= 0; --k) {
                        if( taxonSets.get(k).contains(c) ) {
                            if( ! candidateClades.contains(k) ) {
                                if( onParent.get(k) ) {
                                    if( !intersecting(m_bounds.get(k), m_bounds.get(i)) ) {
                                        break;
                                    }
                                } else {
                                  if(m_bounds.get(k).lower - 1e10 > m_bounds.get(i).lower) {
                                      break;
                                  }
                                }
                                candidateClades.add(k);
                            }
                            break;
                        }
                        if( k == 0 ) {
                           canTaxa.add(c);
                        }
                    }
                }

                final int sz1 = canTaxa.size();
                final int sz2 = candidateClades.size();

                if( sz1 + sz2 == 0 && i + 1 == enclosingIndex  ) {
                    final Bound ebound = m_bounds.get(enclosingIndex);
                    ebound.restrict(m_bounds.get(i));
                } else {
                	if (sz1 + sz2 == 0) {
                		int h = 3;
                				h++;
                	}
                    assert sz1 + sz2 > 0;
                    // prefer taxa over clades (less chance of clades useOriginate clashing)
                    final int k = Randomizer.nextInt(sz1 > 0 ? sz1 : sz2);
                    Set<String> connectTo;
                    int insertPoint;
                    if( k < sz1 ) {
                        // from taxa
                        connectTo = new HashSet<>(1);
                        connectTo.add(canTaxa.get(k));
                        insertPoint = i + 1;
                    } else {
                        // from clade
                        final Iterator<Integer> it = candidateClades.iterator();
                        for (j = 0; j < k - sz1 - 1; ++j) {
                            it.next();
                        }
                        insertPoint = it.next();
                        connectTo = new HashSet<>(taxonSets.get(insertPoint));
                        insertPoint = Math.max(insertPoint, i) + 1;
                    }

                    final HashSet<String> cc = new HashSet<String>(connectTo);

                    connectTo.addAll(taxonSets.get(i));
                    if( !connectTo.equals(enclosingTaxa) || enclosingTaxa == sTaxa ) { // equal when clade already exists

                        taxonSets.add(insertPoint, connectTo);
                        distributions.add(insertPoint, distributions.get(i));
                        onParent.add(insertPoint, false);
                        m_bounds.add(insertPoint, m_bounds.get(i));
                        final String tid = taxonSetIDs.get(i);
                        taxonSetIDs.add(insertPoint, tid);
                        lastMonophyletic += 1;
                    } else {
                        // we lose distribution i :(
                        final Bound ebound = m_bounds.get(enclosingIndex);
                        ebound.restrict(m_bounds.get(i));
                    }
                }
                if( true ) {
                    taxonSets.set(i, new HashSet<>());
                    distributions.set(i, null);
                    m_bounds.set(i, new Bound());
                    final String tid = taxonSetIDs.get(i);
                    if( tid != null ) {
                        taxonSetIDs.set(i, "was-" + tid);
                    }
                }
            }
        }

        {
            int icur = 0;
            for (int i = 0; i < lastMonophyletic; ++i, ++icur) {
                final Set<String> ti = taxonSets.get(i);
                if( ti.isEmpty() ) {
                    icur -= 1;
                } else {
                    if( icur < i ) {
                        taxonSets.set(icur, taxonSets.get(i));
                        distributions.set(icur, distributions.get(i));
                        m_bounds.set(icur, m_bounds.get(i));
                        taxonSetIDs.set(icur, taxonSetIDs.get(i));
                        onParent.set(icur, onParent.get(i));
                    }
                }
            }
            taxonSets.subList(icur, lastMonophyletic).clear();
            distributions.subList(icur, lastMonophyletic).clear();
            m_bounds.subList(icur, lastMonophyletic).clear();
            taxonSetIDs.subList(icur, lastMonophyletic).clear();
            onParent.subList(icur, lastMonophyletic).clear();

            lastMonophyletic = icur;
        }

        if( ICC ) {
            for (int i = 0; i < lastMonophyletic; i++) {
                final Set<String> ti = taxonSets.get(i);
                for (int j = i + 1; j < lastMonophyletic; j++) {
                    final Set<String> tj = taxonSets.get(j);
                    boolean ok = tj.containsAll(ti);
                    if( ok ) {
                        ok = !tj.equals(ti) || (!onParent.get(i) && onParent.get(j));
                        assert ok : "" + i + ' ' + j + ' ' + taxonSetIDs.get(i) + ' ' + taxonSetIDs.get(j);
                    } else {
                        Set<String> tmp = new HashSet<>(tj);
                        tmp.retainAll(ti);
                        assert tmp.isEmpty();
                    }
                }
            }
        }

        // map parent child relationships between mono clades. nParent[i] is the immediate parent clade of i, if any. An immediate parent is the
        // smallest superset of i, children[i] is a list of all clades which have i as a parent.
        // The last one, standing for the virtual "root" of all monophyletic clades is not associated with any actual clade
        final int[] nParent = new int[lastMonophyletic];
        children = new List[lastMonophyletic + 1];
        for (int i = 0; i < lastMonophyletic + 1; i++) {
            children[i] = new ArrayList<>();
        }
        for (int i = 0; i < lastMonophyletic; i++) {
            int j = i + 1;
            while (j < lastMonophyletic && !taxonSets.get(j).containsAll(taxonSets.get(i))) {
                j++;
            }
            nParent[i] = j;
            children[j].add(i);
        }

        // make sure upper bounds of a child does not exceed the upper bound of its parent
        for (int i = lastMonophyletic-1; i >= 0 ;--i) {
            if (nParent[i] < lastMonophyletic ) {
                if (m_bounds.get(i).upper > m_bounds.get(nParent[i]).upper) {
                    m_bounds.get(i).upper = m_bounds.get(nParent[i]).upper - 1e-100;
                    
                    if (m_bounds.get(i).lower >  m_bounds.get(i).upper) {
                    	int h = 3;
                    	h++;
                    }
                    assert m_bounds.get(i).lower <=  m_bounds.get(i).upper: i;
                }
            }
        }

        nodeCount = 2 * sTaxa.size() - 1;
        boundPerNode = new Bound[nodeCount];
        distPerNode = new ParametricDistribution[nodeCount];

        buildTree(sTaxa);                                         
        assert nextNodeNr == nodeCount : "" + nextNodeNr + ' ' + nodeCount;

        double bm = branchMeanInput.get();

        if( bm < 0 ) {
            double maxMean = 0;

            for (ParametricDistribution distr : distPerNode) {
                if( distr != null ) {
                    double m = distr.getMean();
                    if( maxMean < m ) maxMean = m;
                }
            }
            if( maxMean > 0 ) {
                double s = 0;
                for (int i = 2; i <= nodeCount; ++i) {
                    s += 1.0 / i;
                }
                bm = s / maxMean;
            }
        }

        double rate = 1 / (bm < 0 ? 1 : bm);
        boolean succ = false;
        int ntries = 6;
        final double epsi = 0.01/rate;
        double clamp = 1-clampInput.get();
        while( !succ && ntries > 0 ) {
            try {
				succ = setHeights(rate, false, epsi, clamp);
			} catch (ConstraintViolatedException e) {
				throw new RuntimeException("Constraint failed setting heights: " + e.getMessage() + " Perhaps setting/reducing 'branchMean' helps");
			}
            --ntries;
            rate *= 2;
            clamp /= 2;
        }
        if( ! succ ) {
           try {
        	   succ = setHeights(rate, true, 0, 0);
           } catch (ConstraintViolatedException e) {
				throw new RuntimeException("Constraint failed setting heights: " + e.getMessage() + " Perhaps setting/reducing 'branchMean' helps");
           }
        }
        assert succ;

        internalNodeCount = sTaxa.size() - 1;
        leafNodeCount = sTaxa.size();

        HashMap<String,Integer> taxonToNR = null;
        // preserve node numbers where possible
        if (m_initial.get() != null) {
            taxonToNR = new HashMap<>();
            for( Node n : m_initial.get().getExternalNodes() ) {
                taxonToNR.put(n.getID(), n.getNr());
            }
        }
        // re-assign node numbers
        setNodesNrs(root, 0, new int[1], taxonToNR);

        initArrays();
    }


	protected boolean setHeights(final double rate, final boolean safe, final double epsi,
			final double clampBoundsLevel) throws ConstraintViolatedException {
		// node low >= all child nodes low. node high < parent high
		assert rate > 0;
		assert 0 <= clampBoundsLevel && clampBoundsLevel < 1;

		Node[] post = listNodesPostOrder(null, null);
		postCache = null; // can't figure out the javanese to call
							// TreeInterface.listNodesPostOrder

		Bound[] bounds = new Bound[boundPerNode.length];

		for (int i = post.length - 1; i >= 0; --i) {
			final Node node = post[i];
			final int nr = node.getNr();
			bounds[nr] = boundPerNode[nr];

			Bound b = bounds[nr];
			if (b == null) {
				bounds[nr] = b = new Bound();
			}
			if (b.lower < 0) {
				b.lower = 0.0;
			}
			if (clampBoundsLevel > 0) {
				final ParametricDistribution distr = distPerNode[nr];
				if (distr != null) {
					try {
						final double low = distr.inverseCumulativeProbability(clampBoundsLevel / 2);
						final double high = distr.inverseCumulativeProbability(1 - clampBoundsLevel / 2);
						if (distr.density(low) != distr.density(distr.getMean())) {
							if (b.upper >= low && high >= b.lower) {
								b.lower = Math.max(b.lower, low);
								b.upper = Math.min(b.upper, high);
							}
						}
					} catch (MathException e) {
					} catch (RuntimeException e) {
						b.lower = 0.0;
						b.upper = Double.POSITIVE_INFINITY;
						// e.printStackTrace();
					}
				}
			}

			final Node p = node.getParent();
			if (p != null) {
				// assert p.getNr() < bounds.length;
				Bound pb = bounds[p.getNr()];
				assert (pb != null);
				if (!pb.upper.isInfinite() && !(b.upper < pb.upper)) {
					b.upper = pb.upper;
				}
			}
		}

		for (Node node : post) {
			final int nr = node.getNr();
			if (node.isLeaf()) {
				// Bound b = boundPerNode[nr]; assert b != null;
			} else {
				Bound b = bounds[nr];
				assert (b != null);
				for (Node c : node.getChildren()) {
					final Bound cbnd = bounds[c.getNr()];
					b.lower = Math.max(b.lower, cbnd.lower);
					cbnd.upper = Math.min(cbnd.upper, b.upper);
				}
				if (b.lower > b.upper) {
					//throw new ConstraintViolatedException(b.lower + " >" + b.upper);
				}
			}
		}

		if (rootHeightInput.get() != null) {
			final double h = rootHeightInput.get();
			Bound b = bounds[root.getNr()];
			if (b.lower <= h && h <= b.upper) {
				b.upper = h;
			}
		}

		for (Node node : post) {
			if (!node.isLeaf()) {
				final int nr = node.getNr();
				Bound b = bounds[nr];
				double h = -1;
				for (Node c : node.getChildren()) {
					h = Math.max(c.getHeight(), h);
				}
				if (h > b.upper) {
					//throw new ConstraintViolatedException(h + " > " + b);
				}
				if (b.upper.isInfinite()) {
					if (!b.lower.isInfinite()) {
						h = Math.max(h, b.lower);
					}
					h += Randomizer.nextExponential(rate);
				} else {
					if (!b.lower.isInfinite()) {
						h = Math.max(b.lower, h);
					}

					final double range = b.upper - h;
					double r;
					if (safe) {
						r = (range / post.length) * Randomizer.nextDouble();
						assert r > 0 && h + r < b.upper;
					} else {
						r = Randomizer.nextExponential(rate);
						if (r >= range) {
							r = range * Randomizer.nextDouble();
						}
					}
					assert h + r <= b.upper;
					if (r <= epsi && h + r + epsi * 1.001 < b.upper) {
						r += 1.001 * epsi;
					}
					h += r;
				}
				node.setHeight(h);
			}
		}

		if (rootHeightInput.get() != null) {
			final double h = rootHeightInput.get();
			root.setHeight(h);
		}

		// for now fail - this happens rarely
		for (int i = post.length - 1; i >= 0; --i) {
			final Node node = post[i];
			if (!node.isRoot() && node.getLength() <= epsi) {
				return false;
			}
		}
		return true;
	}
}
