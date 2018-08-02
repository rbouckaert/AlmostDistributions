package beast.math.distributions;

import org.apache.commons.math.MathException;

import beast.core.Input;

public class AlmostLogNormal extends LogNormalDistributionModel implements AlmostDistribution {
	public Input<Double> penaltyInput = new Input<Double>("penalty", "penalty for being outside range", 10000.0);

	double _lower, _upper, penalty;
	
	
	@Override
	public void initAndValidate() {
		dist = new AlmostLogNormalImpl(0.0, 1.0);
		
		super.initAndValidate();
		
		penalty = penaltyInput.get();
		
		try {
			_lower = dist.inverseCumulativeProbability(0.025);
			_upper = dist.inverseCumulativeProbability(0.975);
		} catch (MathException e) {
			e.printStackTrace();
			throw new IllegalArgumentException(e);
		}
	}
	
	class AlmostLogNormalImpl extends LogNormalImpl {
        public AlmostLogNormalImpl(double mean, double stdDev) {
			super(mean, stdDev);
		}

		@Override
        public double logDensity(final double x) {
            if (x >= _lower && x <= _upper) {
            	return Math.log(density(x));
            } else if (x < _lower) {
                return -penalty * (_lower - x);
            }  else {
                return -penalty * (x - _upper);
            }
        }
    } // class AlmostLogNormalImpl
	
	@Override
	public double getLowerTarget() throws MathException {
		return _lower;
	}
	
	@Override
	public double getUpperTarget() throws MathException {
		return _upper;
	}

}
