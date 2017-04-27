package com.graf.model;



import java.util.HashMap;
import java.util.Map;

import Modelandprocess.AbstractModelBates;
import Modelandprocess.AbstractProcessInterfaceBates;
import net.finmath.montecarlo.BrownianMotion;
import net.finmath.montecarlo.RandomVariable;
import net.finmath.montecarlo.assetderivativevaluation.AssetModelMonteCarloSimulationInterface;
import net.finmath.stochastic.RandomVariableInterface;
import net.finmath.time.TimeDiscretizationInterface;
import net.finmath.exception.CalculationException;

public class MonteCarloTwoFactorBatesModel extends AbstractModelBates implements AssetModelMonteCarloSimulationInterface {

	/**
	 * Truncation schemes to be used in the calculation of drift and diffusion coefficients.
	 */
	public enum Scheme {
		/**
		 * Reflection scheme, that is V is replaced by Math.abs(V), where V denotes the current realization of V(t).
		 */
		REFLECTION,
		
		/**
		 * Full truncation scheme, that is V is replaced by Math.max(V,0), where V denotes the current realization of V(t).
		 */
		FULL_TRUNCATION
	};

	private final int numberOfPaths;
	private final int seed;
	private final double initialValue;
	private final double riskFreeRate;		// Actually the same as the drift (which is not stochastic) /mu
	private final double volatilityForVOne;
	private final double volatilityForVTwo;
	private final double cvOne;
	private final double cvTwo;

	private final double alphaOne;
	private final double betaOne;
	private final double alphaTwo;
	private final double betaTwo;
	private final double sigmaOne;
	private final double sigmaTwo;
	private final double rhoOne;
	private final double rhoTwo;
	
	
	private final double lambdaZero;
	private final double lambdaOne;
	private final double lambdaTwo;
	private final double k;
	private final double delta;

	private final Scheme scheme;
	
	

	/*
	 * The interface definition requires that we provide the initial value, the drift and the volatility in terms of random variables.
	 * We construct the corresponding random variables here and will return (immutable) references to them.
	 */
	private RandomVariableInterface[]	initialValueVector	= new RandomVariableInterface[3];
	private TimeDiscretizationInterface timeDiscretization;
	
	public MonteCarloTwoFactorBatesModel(
			TimeDiscretizationInterface timeDiscretization,
			int numberOfPaths,
			int seed,
			double initialValue,
			double riskFreeRate,
			double volatilityForVOne,
			double volatilityForVTwo,
			double cvOne,
			double cvTwo,
			double alphaOne,
			double betaOne,
			double alphaTwo,
			double betaTwo,
			double sigmaOne,
			double sigmaTwo,
			double rhoOne,
			double rhoTwo,
			double lambdaZero,
			double lambdaOne,
			double lambdaTwo,
			double k,
			double delta,
			Scheme scheme
			) {
		super();

		this.timeDiscretization = timeDiscretization;
		this.numberOfPaths		= numberOfPaths;
		this.seed				= seed;
		this.initialValue		= initialValue;
		this.riskFreeRate		= riskFreeRate;
		this.volatilityForVOne	= volatilityForVOne;
		this.volatilityForVTwo	= volatilityForVTwo;
		this.cvOne				= cvOne;
		this.cvTwo				= cvTwo;
		this.alphaOne			= alphaOne;
		this.betaOne			= betaOne;
		this.alphaTwo			= alphaTwo;
		this.betaTwo			= betaTwo;
		this.sigmaOne			= sigmaOne;
		this.sigmaTwo			= sigmaTwo;
		this.rhoOne				= rhoOne;
		this.rhoTwo				= rhoTwo;
		this.lambdaZero			= lambdaZero;
		this.lambdaOne			= lambdaOne;
		this.lambdaTwo			= lambdaTwo;
		this.k					= k;
		this.delta				= delta;
		this.scheme				= scheme;
		this.timeDiscretization = timeDiscretization;
		


		BrownianMotion uncorrelatedFactors = new BrownianMotion(timeDiscretization, 4 , numberOfPaths, seed);
		

		// Create a corresponding MC process
		AbstractProcessInterfaceBates process = new ProcessEulerSchemeBates(uncorrelatedFactors);

		// Link model and process for delegation
		process.setModel(this);
		this.setProcess(process);
	}




	public RandomVariableInterface[] getInitialState() {
		// Since the underlying process is configured to simulate log(S), the initial value and the drift are transformed accordingly.
		if(initialValueVector[0] == null) 	{
			initialValueVector[0] = getRandomVariableForConstant(Math.log(initialValue));
			initialValueVector[1] = getRandomVariableForConstant(volatilityForVOne*volatilityForVOne);
			initialValueVector[2] = getRandomVariableForConstant(volatilityForVTwo*volatilityForVTwo);
		}

		return initialValueVector;
	}

	public RandomVariableInterface[] getDrift(int timeIndex, RandomVariableInterface[] realizationAtTimeIndex, RandomVariableInterface[] realizationPredictor) {
		RandomVariableInterface stochasticVarianceOne = realizationAtTimeIndex[1];
//		if(scheme == Scheme.FULL_TRUNCATION)	stochasticVarianceOne = realizationAtTimeIndex[1].floor(0.0);
//		else if(scheme == Scheme.REFLECTION)	stochasticVarianceOne = realizationAtTimeIndex[1].abs();
//		else throw new UnsupportedOperationException("Scheme " + scheme.name() + " not supported.");
		
		RandomVariableInterface stochasticVarianceTwo = realizationAtTimeIndex[2];
//		if(scheme == Scheme.FULL_TRUNCATION)	stochasticVarianceTwo = realizationAtTimeIndex[2].floor(0.0);
//		else if(scheme == Scheme.REFLECTION)	stochasticVarianceTwo = realizationAtTimeIndex[2].abs();
//		else throw new UnsupportedOperationException("Scheme " + scheme.name() + " not supported.");
		
		//calculate time-varying lambda
		RandomVariableInterface lambda = getRandomVariableForConstant(lambdaZero).addProduct(stochasticVarianceOne, lambdaOne).addProduct(stochasticVarianceTwo, lambdaTwo);

		RandomVariableInterface[] drift = new RandomVariable[3];
		
		//-1/2 (V1+V2) because of state space transform for the drift[0]
		drift[0] = getRandomVariableForConstant(riskFreeRate).addProduct(lambda, -k).addProduct(stochasticVarianceOne, cvOne-0.5).addProduct(stochasticVarianceTwo, cvTwo-0.5);
		drift[1] = getRandomVariableForConstant(alphaOne).addProduct(stochasticVarianceOne, -betaOne);
		drift[2] = getRandomVariableForConstant(alphaTwo).addProduct(stochasticVarianceTwo, -betaTwo);
		return drift;
	}

	public RandomVariableInterface[] getFactorLoading(int timeIndex, int component, RandomVariableInterface[] realizationAtTimeIndex) {
		RandomVariableInterface stochasticVarianceOne;
		if(scheme == Scheme.FULL_TRUNCATION)	stochasticVarianceOne = realizationAtTimeIndex[1].floor(0.0).sqrt();
		else if(scheme == Scheme.REFLECTION)	stochasticVarianceOne = realizationAtTimeIndex[1].abs().sqrt();
		else throw new UnsupportedOperationException("Scheme " + scheme.name() + " not supported.");
		
		RandomVariableInterface stochasticVarianceTwo;
		if(scheme == Scheme.FULL_TRUNCATION)	stochasticVarianceTwo = realizationAtTimeIndex[2].floor(0.0).sqrt();
		else if(scheme == Scheme.REFLECTION)	stochasticVarianceTwo = realizationAtTimeIndex[2].abs().sqrt();
		else throw new UnsupportedOperationException("Scheme " + scheme.name() + " not supported.");

		RandomVariableInterface[] factorLoadings = new RandomVariableInterface[4];

		if(component == 0) {
			factorLoadings[0] = stochasticVarianceOne;
			factorLoadings[1] = getRandomVariableForConstant(0.0);
			factorLoadings[2] = stochasticVarianceTwo;
			factorLoadings[3] = getRandomVariableForConstant(0.0);
		}
		else if(component == 1) {
			factorLoadings[0] = stochasticVarianceOne.mult(sigmaOne).mult(rhoOne);
			factorLoadings[1] = stochasticVarianceOne.mult(sigmaOne).mult(Math.sqrt(1-rhoOne*rhoOne));
			factorLoadings[2] = getRandomVariableForConstant(0.0);
			factorLoadings[3] = getRandomVariableForConstant(0.0);
		}
		else if(component == 2) {
			factorLoadings[0] = getRandomVariableForConstant(0.0);
			factorLoadings[1] = getRandomVariableForConstant(0.0);
			factorLoadings[2] = stochasticVarianceTwo.mult(sigmaTwo).mult(rhoTwo);
			factorLoadings[3] = stochasticVarianceTwo.mult(sigmaTwo).mult(Math.sqrt(1-rhoTwo*rhoTwo));
		}
		else {
			throw new UnsupportedOperationException("Component " + component + " does not exist.");
		}

		return factorLoadings;
	}
	
	@Override
	public RandomVariableInterface getJumpComponent(int timeIndex, RandomVariableInterface[] realizationAtTimeIndex) {
		RandomVariableInterface stochasticVarianceOne;
		if(scheme == Scheme.FULL_TRUNCATION)	stochasticVarianceOne = realizationAtTimeIndex[1].floor(0.0);
		else if(scheme == Scheme.REFLECTION)	stochasticVarianceOne = realizationAtTimeIndex[1].abs();
		else throw new UnsupportedOperationException("Scheme " + scheme.name() + " not supported.");
		
		RandomVariableInterface stochasticVarianceTwo;
		if(scheme == Scheme.FULL_TRUNCATION)	stochasticVarianceTwo = realizationAtTimeIndex[2].floor(0.0);
		else if(scheme == Scheme.REFLECTION)	stochasticVarianceTwo = realizationAtTimeIndex[2].abs();
		else throw new UnsupportedOperationException("Scheme " + scheme.name() + " not supported.");
		
		//calculate time-varying lambda
		RandomVariableInterface lambda = getRandomVariableForConstant(lambdaZero).addProduct(stochasticVarianceOne, lambdaOne).addProduct(stochasticVarianceTwo, lambdaTwo);
		
		
		
		return PoissonBates.getPoisson(timeDiscretization.getTimeStep(timeIndex), numberOfPaths, k, delta, lambda, seed);
	}

	public RandomVariableInterface applyStateSpaceTransform(int componentIndex, RandomVariableInterface randomVariable) {
		if(componentIndex == 0) {
			return randomVariable.exp();
		}
		else if(componentIndex == 1) {
			return randomVariable;
		}
		else if(componentIndex == 2) {
			return randomVariable;
		}
		else {
			throw new UnsupportedOperationException("Component " + componentIndex + " does not exist.");
		}
	}

	public int getNumberOfComponents() {
		return 3;
	}

	public RandomVariableInterface getRandomVariableForConstant(double value) {
		return getProcess().getStochasticDriver().getRandomVariableForConstant(value);
	}

	public MonteCarloTwoFactorBatesModel getCloneWithModifiedData(Map<String, Object> dataModified) {
		/*
		 * Determine the new model parameters from the provided parameter map.
		 */
		double	newInitialTime			= dataModified.get("initialTime") != null		? ((Number)dataModified.get("initialTime")).doubleValue() : getTime(0);
		int		newSeed					= dataModified.get("seed") != null				? ((Number)dataModified.get("seed")).intValue() : seed;
		double	newInitialValue			= dataModified.get("initialValue") != null		? ((Number)dataModified.get("initialValue")).doubleValue() : initialValue;
		double	newRiskFreeRate			= dataModified.get("riskFreeRate") != null		? ((Number)dataModified.get("riskFreeRate")).doubleValue() : riskFreeRate;
		double	newVolatilityForVOne	= dataModified.get("volatilityForVOne") != null	? ((Number)dataModified.get("volatilityForVOne")).doubleValue() : volatilityForVOne;
		double	newVolatilityForVTwo	= dataModified.get("volatilityForVTwo") != null	? ((Number)dataModified.get("volatilityForVTwo")).doubleValue() : volatilityForVTwo;
		double	newCvOne				= dataModified.get("cvOne") != null				? ((Number)dataModified.get("cvOne")).doubleValue() : cvOne;
		double	newCvTwo				= dataModified.get("cvTwo") != null				? ((Number)dataModified.get("cvTwo")).doubleValue() : cvTwo;
		double	newAlphaOne				= dataModified.get("alphaOne") != null			? ((Number)dataModified.get("alphaOne")).doubleValue() : alphaOne;
		double	newBetaOne				= dataModified.get("betaOne") != null			? ((Number)dataModified.get("betaOne")).doubleValue() : betaOne;
		double	newAlphaTwo				= dataModified.get("alphaTwo") != null			? ((Number)dataModified.get("alphaTwo")).doubleValue() : alphaTwo;
		double	newBetaTwo				= dataModified.get("betaTwo") != null			? ((Number)dataModified.get("betaTwo")).doubleValue() : betaTwo;
		double	newSigmaOne				= dataModified.get("sigmaOne") != null			? ((Number)dataModified.get("sigmaOne")).doubleValue() : sigmaOne;
		double	newSigmaTwo				= dataModified.get("sigmaTwo") != null			? ((Number)dataModified.get("sigmaTwo")).doubleValue() : sigmaTwo;
		double	newRhoOne				= dataModified.get("rhoOne") != null			? ((Number)dataModified.get("rhoOne")).doubleValue() : rhoOne;
		double	newRhoTwo				= dataModified.get("rhoTwo") != null			? ((Number)dataModified.get("rhoTwo")).doubleValue() : rhoTwo;
		double	newLambdaZero			= dataModified.get("lambdaZero") != null		? ((Number)dataModified.get("lambdaZero")).doubleValue() : lambdaZero;
		double	newLambdaOne			= dataModified.get("lambdaOne") != null			? ((Number)dataModified.get("lambdaOne")).doubleValue() : lambdaOne;
		double	newLambdaTwo			= dataModified.get("lambdaTwo") != null			? ((Number)dataModified.get("lambdaTwo")).doubleValue() : lambdaTwo;
		double	newK					= dataModified.get("k") != null					? ((Number)dataModified.get("k")).doubleValue() : k;
		double	newDelta				= dataModified.get("delta") != null				? ((Number)dataModified.get("delta")).doubleValue() : delta;
		
		return new MonteCarloTwoFactorBatesModel(
				getProcess().getTimeDiscretization().getTimeShiftedTimeDiscretization(newInitialTime-getTime(0)),
				getProcess().getNumberOfPaths(),
				newSeed,
				newInitialValue,
				newRiskFreeRate,
				newVolatilityForVOne,
				newVolatilityForVTwo,
				newCvOne,
				newCvTwo,
				newAlphaOne,
				newBetaOne,
				newAlphaTwo,
				newBetaTwo,
				newSigmaOne,
				newSigmaTwo,
				newRhoOne,
				newRhoTwo,
				newLambdaZero,
				newLambdaOne,
				newLambdaTwo,
				newK,
				newDelta,
				scheme
				);
	}

	@Override
	public String toString() {
		return "TwoFactorBatesModel [initialValue=" + initialValue + ", riskFreeRate=" + riskFreeRate + ", volatilityForVOne="
				+ volatilityForVOne + ", volatilityForVTwo=" + volatilityForVTwo + ", cvOne=" + cvOne + ", cvTwo=" + cvTwo + ", alphaOne=" + alphaOne 
				+ ", betaOne=" + betaOne + ", alphaTwo=" + alphaTwo + ", betaTwo=" + betaTwo 
				+ ", sigmaOne=" + sigmaOne + ", sigmaTwo=" + sigmaTwo + ", rhoOne=" + rhoOne
				+ ", rhoTwo=" + rhoTwo + ", lambdaZero=" + lambdaZero + ", lambdaOne=" + lambdaOne
				+ ", lambdaTwo=" + lambdaTwo + ", k=" + k + ", delta=" + delta + ", scheme="
				+ scheme + "]";
	}

	/**
	 * Returns the risk free rate parameter of this model.
	 *
	 * @return Returns the riskFreeRate.
	 */
	public double getRiskFreeRate() {
		return riskFreeRate;
	}

	/**
	 * Returns the volatility parameter of this model.
	 * 
	 * @return Returns the volatility.
	 */
//	public double getVolatility() {
//		return volatility;
//	}
	@Override
	public RandomVariableInterface getAssetValue(double time, int assetIndex) throws CalculationException {
		return getAssetValue(getTimeIndex(time), assetIndex);
	}

	@Override
	public RandomVariableInterface getAssetValue(int timeIndex, int assetIndex) throws CalculationException {
		return getProcess().getProcessValue(timeIndex, assetIndex);
	}



	@Override
	public RandomVariableInterface getMonteCarloWeights(double time) throws CalculationException {
		return getMonteCarloWeights(getTimeIndex(time));
	}

	/* (non-Javadoc)
	 * @see net.finmath.montecarlo.assetderivativevaluation.AssetModelMonteCarloSimulationInterface#getNumberOfAssets()
	 */
	@Override
	public int getNumberOfAssets() {
		return 1;
	}

	@Override
	public int getNumberOfPaths() {

		return numberOfPaths;
	}

	@Override
	public RandomVariableInterface getNumeraire(int timeIndex) throws CalculationException {
		double time = getTime(timeIndex);

		return getNumeraire(time);
	}

	@Override
	public RandomVariableInterface getNumeraire(double time) throws CalculationException {
		return getProcess().getStochasticDriver().getRandomVariableForConstant(Math.exp(riskFreeRate * time));
	}

	@Override
	public AssetModelMonteCarloSimulationInterface getCloneWithModifiedSeed(int seed) throws CalculationException {
		Map<String, Object> dataModified = new HashMap<String, Object>();
		dataModified.put("seed", new Integer(seed));
		return getCloneWithModifiedData(dataModified);
	}



}

	
	

