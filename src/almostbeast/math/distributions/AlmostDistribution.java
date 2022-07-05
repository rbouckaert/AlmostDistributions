package almostbeast.math.distributions;

import org.apache.commons.math.MathException;

import beast.base.inference.distribution.ParametricDistribution;

public interface AlmostDistribution {

	default public double getLowerTarget() throws MathException {
		return ((ParametricDistribution)this).inverseCumulativeProbability(0.05);
	}

	default public double getUpperTarget() throws MathException {
		return ((ParametricDistribution)this).inverseCumulativeProbability(0.95);
	}
}
