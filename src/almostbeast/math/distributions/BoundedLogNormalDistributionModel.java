package almostbeast.math.distributions;


import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.ContinuousDistribution;
import org.apache.commons.math.distribution.Distribution;
import org.apache.commons.math.distribution.IntegerDistribution;

import beast.base.core.Input;
import beast.base.inference.distribution.LogNormalDistributionModel;

public class BoundedLogNormalDistributionModel extends LogNormalDistributionModel {
	public Input<Double> lowerInput = new Input<Double>("lower","lower bound on the interval, defaul 0", 0.0);
	public Input<Double> upperInput = new Input<Double>("upper","lower bound on the interval, defaul 1", 1.0);

	double lower, upper;
	
	@Override
	public void initAndValidate() {
		super.initAndValidate();
		lower = lowerInput.get();
		upper = upperInput.get();
		if (lower > upper) {
			throw new IllegalArgumentException("Lower should be higher than upper");
		}
	}

	
    public double logDensity(double x) {
        x -= offsetInput.get();
        if (x < lower || x > upper) {
            return -10000 * (x-getMean()) * (x-getMean());
        	// return Double.NEGATIVE_INFINITY;
        }
        final Distribution dist = getDistribution();
        if (dist instanceof ContinuousDistribution) {
            return ((ContinuousDistribution) dist).logDensity(x);
        } else if (dist instanceof IntegerDistribution) {
            final double probability = ((IntegerDistribution) dist).probability(x);
            if( probability > 0 ) {
                return Math.log(probability);
            }
        }
        return Double.NEGATIVE_INFINITY;
    }

    @Override
	public double cumulativeProbability(final double x) throws MathException {
    	throw new RuntimeException("Not implemented yet");
        // return getDistribution().cumulativeProbability(x);
    }

    @Override
	public double cumulativeProbability(final double x0, final double x1) throws MathException {
    	throw new RuntimeException("Not implemented yet");
        // return getDistribution().cumulativeProbability(x0, x1);
    }

    @Override
	public double inverseCumulativeProbability(final double p) throws MathException {
        double offset = offsetInput.get();
        if (offset != 0) {
        	throw new RuntimeException("Not implemented yet when offset != 0.0");
        }
        final Distribution dist = getDistribution();
        double plo = dist.cumulativeProbability(lower);
        double pup = dist.cumulativeProbability(upper);
        double diff = pup - plo;
        return offset + ((ContinuousDistribution) dist).inverseCumulativeProbability(plo + p * diff);
    }

}
