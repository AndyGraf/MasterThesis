/*
 * (c) Copyright Christian P. Fries, Germany. All rights reserved. Contact: email@christian-fries.de.
 *
 * Created on 23.03.2014
 */
package com.graf.fouriertransform;

import org.apache.commons.math3.complex.Complex;

import net.finmath.compatibility.java.util.function.DoubleUnaryOperator;
import net.finmath.fouriermethod.CharacteristicFunctionInterface;
import net.finmath.fouriermethod.models.ProcessCharacteristicFunctionInterface;
import net.finmath.fouriermethod.products.AbstractProductFourierTransform;
import net.finmath.integration.RealIntegralInterface;
import net.finmath.integration.SimpsonRealIntegrator;

/**
 * Implements valuation of a European option on a single asset.
 * 
 * Given a model for an asset <i>S</i>, the European option with strike <i>K</i>, maturity <i>T</i>
 * pays
 * <br>
 * 	<i>max(S(T) - K , 0)</i> in <i>T</i>
 * <br>
 * 
 * The class implements the characteristic function of the call option
 * payoff, i.e., its Fourier transform.
 * 
 * @author Christian Fries
 * @author Alessandro Gnoatto
 * @version 1.0
 */
public class EuropeanOption extends AbstractProductFourierTransform {

	private final double maturity;
	private final double strike;
	private final double payoffUnit;
	private final String nameOfUnderliyng;
	
	/**
	 * Construct a product representing an European option on an asset S (where S the asset with index 0 from the model - single asset case).
	 * @param maturity The maturity T in the option payoff max(S(T)-K,0)
	 * @param strike The strike K in the option payoff max(S(T)-K,0).
	 */
	public EuropeanOption(double maturity, double strike, double payoffUnit) {
		super();
		this.maturity			= maturity;
		this.strike				= strike;
		this.payoffUnit			= payoffUnit;
		this.nameOfUnderliyng	= null;		// Use asset with index 0
	}

	/* (non-Javadoc)
	 * @see net.finmath.fouriermethod.CharacteristicFunctionInterface#apply(org.apache.commons.math3.complex.Complex)
	 */
	@Override
	public Complex apply(Complex argument) {
		Complex iargument = argument.multiply(Complex.I);
		Complex exponent = (iargument).add(1);
		Complex numerator = (new Complex(strike)).pow(exponent);
		Complex denominator = (argument.multiply(argument)).subtract(iargument);	

		return numerator.divide(denominator).negate();
	}


	//includes a discounting factor
	@Override
	public double getValue(ProcessCharacteristicFunctionInterface model) {

		final CharacteristicFunctionInterface modelCF   = model.apply(getMaturity());
        final CharacteristicFunctionInterface productCF = this;

		final double lineOfIntegration = 0.5 * getIntegrationDomainImagUpperBound()+getIntegrationDomainImagLowerBound();
		DoubleUnaryOperator integrand = new DoubleUnaryOperator() {
            @Override
            public double applyAsDouble(double real) {
                Complex z = new Complex(real,lineOfIntegration);
                return modelCF.apply(z.negate()).multiply(productCF.apply(z)).getReal();
            }
        };

		RealIntegralInterface integrator = new SimpsonRealIntegrator(-100.0, 100.0, 20000, true);

		return payoffUnit * integrator.integrate(integrand) / 2.0 / Math.PI;
	}
	
	
	/* (non-Javadoc)
	 * @see net.finmath.fouriermethod.products.AbstractProductFourierTransform#getMaturity()
	 */
	@Override
	public double getMaturity() {
		return maturity;
	}

	/* (non-Javadoc)
	 * @see net.finmath.fouriermethod.products.AbstractProductFourierTransform#getDomainImagLowerBound()
	 */
	@Override
	public double getIntegrationDomainImagLowerBound() {
		return 0.5;
	}

	/* (non-Javadoc)
	 * @see net.finmath.fouriermethod.products.AbstractProductFourierTransform#getDomainImagUpperBound()
	 */
	@Override
	public double getIntegrationDomainImagUpperBound() {
		return 2.5;
	}
}
