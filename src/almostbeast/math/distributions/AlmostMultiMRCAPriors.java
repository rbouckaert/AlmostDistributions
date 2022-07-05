package almostbeast.math.distributions;


import beast.base.core.Description;
import beast.base.evolution.tree.MRCAPrior;
import beastlabs.math.distributions.MultiMRCAPriors;

@Description("Allows multiple MRCAs to be violated if set to monophyletic")
public class AlmostMultiMRCAPriors extends MultiMRCAPriors {

    @Override
    public double calculateLogP() {
    	logP = super.calculateLogP();
    	if (Double.isInfinite(logP)) {
    		logP = 0;
    		for (MRCAPrior p : calibrationsInput.get()) {
    			// may become finite, if AlmostMRCAPriors are used
    			logP += p.calculateLogP();
    		}
    	}
    	return logP;
    }
}
