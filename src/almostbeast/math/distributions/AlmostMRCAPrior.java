package almostbeast.math.distributions;

import beast.base.core.Description;
import beast.base.evolution.tree.MRCAPrior;

@Description("Allows MRCAs to be violated if set to monophyletic")
public class AlmostMRCAPrior extends MRCAPrior {

	
	@Override
	public double calculateLogP() {
		logP = super.calculateLogP();
		if (Double.isInfinite(logP)) {
			logP = -10000;
		}
		return logP;
	}
}
