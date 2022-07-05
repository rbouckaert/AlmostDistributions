package almostbeast.math.distributions;


import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.Distribution;

import beast.base.inference.distribution.Uniform;
import beast.base.core.Description;
import beast.base.core.Input;

@Description("Uniform distribution with very small support outside range (to prevent infinities)")
public class AlmostUniform extends Uniform implements AlmostDistribution {
	public Input<Double> penaltyInput = new Input<Double>("penalty", "penalty for being outside range", 10000.0);

    double centre;
    double penalty;

    private boolean infiniteSupport;

    @Override
    public void initAndValidate() {
        distr = new AlmostUniformImpl();

        _lower = lowerInput.get();
        _upper = upperInput.get();
        if (_lower >= _upper) {
            throw new IllegalArgumentException("Upper value should be higher than lower value");
        }
        distr.setBounds(_lower, _upper);
        infiniteSupport = Double.isInfinite(_lower) || Double.isInfinite(_upper);
        if (infiniteSupport) {
            density = 1.0;
        } else {
            density = 1.0 / (_upper - _lower);
        }
        
        if (Double.isFinite(_upper) & Double.isFinite(_lower)) {
        	centre = offsetInput.get() + (_upper + _lower)/2;
        } else if (Double.isFinite(_upper)) {
        	centre = _upper;
        } else {
        	centre = _lower;
        }
        
        penalty = penaltyInput.get();
    }


    class AlmostUniformImpl extends UniformImpl {
        private double lower;
        private double upper;

        public void setBounds(final double lower, final double upper) {
            this.lower = lower;
            this.upper = upper;
        }

        @Override
        public double cumulativeProbability(double x) throws MathException {
            x = Math.max(x, lower);
            return (x - lower) / (upper - lower);
        }

        @Override
        public double cumulativeProbability(double x0, double x1) throws MathException {
            x0 = Math.max(x0, lower);
            x1 = Math.min(x1, upper);
            if (x1 < lower || x1 > upper) {
                throw new RuntimeException("Value x (" + x1 + ") out of bounds (" + lower + "," + upper + ").");
            }
            return (x1 - x0) / (upper - lower);
        }

        @Override
        public double inverseCumulativeProbability(final double p) throws MathException {
            if (p < 0.0 || p > 1.0) {
                throw new RuntimeException("inverseCumulativeProbability::argument out of range [0...1]");
            }
            if( p == 0 ) {
                // works even when one bound is infinite
                return _lower;
            }
            if( p == 1 ) {
                // works even when one bound is infinite
                return _upper;
            }
//            if( infiniteSupport ) {
//                 throw new RuntimeException("Inverse Cumulative Probability for 0 < p < 1 and infinite support") ;
//            }
            return (upper - lower) * p + lower;
        }

        @Override
        public double density(final double x) {
            if (x >= lower && x <= upper) {
                return density;
            }
            return 0.0;
        }

        @Override
        public double logDensity(final double x) {
            if (x >= lower && x <= upper) {
            	return Math.log(density(x));
            } else {
                return -penalty * (x-centre) * (x-centre);
            } 
        }
    } // class UniformImpl


    @Override
    public Distribution getDistribution() {
        return distr;
    }

    @Override
    public double density(final double x) {
        if (x >= _lower && x <= _upper) {
            // (BUG)?? why does this not return this.density??? (JH)
            return 1;
        } else {
            return 0;
        }
    }
    
    @Override
    public double getMeanWithoutOffset() {
    	if (Double.isInfinite(_lower) || Double.isInfinite(_upper)) {
    		return Double.NaN;
    	}
    	return offsetInput.get() + (_upper + _lower)/2;
    }

    @Override
    public double getLowerTarget() throws MathException {
    	return _lower;
    }
    
    @Override
    public double getUpperTarget() throws MathException {
    	return _upper;
    }
}
