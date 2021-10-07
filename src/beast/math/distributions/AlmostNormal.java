package beast.math.distributions;


import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.ContinuousDistribution;
import org.apache.commons.math.distribution.NormalDistributionImpl;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;

@Description("Normal distribution with very small support outside range of 2 sigma")
public class AlmostNormal extends ParametricDistribution implements AlmostDistribution {
    final public Input<RealParameter> meanInput = new Input<>("mean", "mean of the normal distribution, defaults to 0");
    final public Input<RealParameter> sigmaInput = new Input<>("sigma", "standard deviation of the normal distribution, defaults to 1");
    final public Input<RealParameter> tauInput = new Input<>("tau", "precission of the normal distribution, defaults to 1", Validate.XOR, sigmaInput);
	final public Input<Double> penaltyInput = new Input<>("penalty", "penalty for being outside range", 10000.0);

	final public Input<Double> upperInput = new Input<>("upper", "upper bound above which a penalty applies", Double.MAX_VALUE);
	final public Input<Boolean> trueBoundsInput = new Input<>("trueBounds", "if true report true upper and lower bounds, if false returns infinite range", true);
	
	AlmostNormalDistributionImpl dist = new AlmostNormalDistributionImpl(0, 1);

    double penalty;
    double customUpper;
    boolean trueBounds = true;

    @Override
    public void initAndValidate() {
    	penalty = penaltyInput.get();
    	customUpper = upperInput.get();
    	trueBounds = trueBoundsInput.get();
    	refresh();
   }
    
	void refresh() {
        double mean;
        double sigma;
        if (meanInput.get() == null) {
            mean = 0;
        } else {
            mean = meanInput.get().getValue();
        }
        if (sigmaInput.get() == null) {
        	if (tauInput.get() == null) {
        		sigma = 1;
        	} else {
                sigma = Math.sqrt(1.0/tauInput.get().getValue());
        	}
        } else {
            sigma = sigmaInput.get().getValue();
        }
        dist.setMean(mean);
        dist.setStandardDeviation(sigma);
    }


    class AlmostNormalDistributionImpl extends NormalDistributionImpl {
 		private static final long serialVersionUID = 1L;

 		private double lower, upper;
    	
    	AlmostNormalDistributionImpl(double mean, double sigma) {
    		super(mean, sigma);
    		setBounds();
    	}
    	
        @SuppressWarnings("deprecation")
        public void setMean(double mean) {
        	super.setMean(mean);
    		setBounds();
        }

        @SuppressWarnings("deprecation")
		public void setStandardDeviation(double sd) {
        	super.setStandardDeviation(sd);
    		setBounds();
        }

        private void setBounds() {
        	lower = this.getMean() - 2 * this.getStandardDeviation();
        	upper = this.getMean() + 2 * this.getStandardDeviation();
        	upper = Math.min(upper, customUpper);
        }

    	
        @Override
        public double cumulativeProbability(double x) throws MathException {
        	if (x < lower) {
        		return 0;
        	}
        	if (x > upper) {
        		return 1;
        	}
        	double p = (super.cumulativeProbability(x) - 0.025)/ 0.95;
        	return p;
        }

        @Override
        public double cumulativeProbability(double x0, double x1) throws MathException {
        	return cumulativeProbability(x1) - cumulativeProbability(x0);
        }

        @Override
        public double inverseCumulativeProbability(final double p) throws MathException {
            if (p < 0.0 || p > 1.0) {
                throw new RuntimeException("inverseCumulativeProbability::argument out of range [0...1]");
            }
            if( p == 0 ) {
            	if (trueBounds)
            		return lower;
                // works even when one bound is infinite
                return Double.NEGATIVE_INFINITY;
            }
            if( p == 1 ) {
            	if (trueBounds)
            		return upper;
                // works even when one bound is infinite
                return Double.POSITIVE_INFINITY;
            }
            double x = super.inverseCumulativeProbability(0.025 + p/ 1.05);
            return x;
        }

        @Override
        public double logDensity(final double x) {
            if (x >= lower && x <= upper) {
            	return Math.log(density(x)/0.95);
            } else {
                return -penalty * (x-getMean()) * (x-getMean());
            } 
        }
    } // class AlmostNormalDistributionImpl

    @Override
    public ContinuousDistribution getDistribution() {
        refresh();
        return dist;
    }

    @Override
    public double getMeanWithoutOffset() {
        if (meanInput.get() == null) {
        	return offsetInput.get();
        } else {
        	return offsetInput.get() + meanInput.get().getValue();
        }
    }
    
    @Override
    public double getLowerTarget() throws MathException {
    	return dist.lower;
    }
    
    @Override
    public double getUpperTarget() throws MathException {
    	return dist.upper;
    }
}
