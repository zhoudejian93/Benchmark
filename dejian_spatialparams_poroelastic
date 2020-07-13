// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup PoromechanicsTests
 * \brief Definition of the spatial parameters for the poro-elastic problem.
 */
#ifndef DUMUX_POROELASTIC_SPATIAL_PARAMS_HH
#define DUMUX_POROELASTIC_SPATIAL_PARAMS_HH

#include <dumux/geomechanics/lameparams.hh>
#include <dumux/material/spatialparams/fvporoelastic.hh>
#include <dumux/material/fluidmatrixinteractions/porositydeformation.hh>

namespace Dumux {

/*!
 * \ingroup PoromechanicsTests
 * \brief Definition of the spatial parameters for the poro-elastic
 *        sub-problem in the coupled poro-mechanical el1p problem.
 */
template<class Scalar, class GridGeometry>
class PoroElasticSpatialParams : public FVSpatialParamsPoroElastic< Scalar,
                                                                    GridGeometry,
                                                                    PoroElasticSpatialParams<Scalar, GridGeometry> >
{
    using ThisType = PoroElasticSpatialParams<Scalar, GridGeometry>;
    using ParentType = FVSpatialParamsPoroElastic<Scalar, GridGeometry, ThisType>;

    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    //! Export the type of the lame parameters
    using LameParams = Dumux::LameParams<Scalar>;

    PoroElasticSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , initPorosity_(getParam<Scalar>("SpatialParams.InitialPorosity"))
    , initUpperAquiferPorosity_(getParam<Scalar>("SpatialParams.InitialUpperAquiferPorosity", 0.15))
    , initCaprockPorosity_(getParam<Scalar>("SpatialParams.InitialCaprockPorosity", 0.01))
    , initReservoirPorosity_(getParam<Scalar>("SpatialParams.InitialReservoirPorosity", 0.2))
    , initBaserockPorosity_(getParam<Scalar>("SpatialParams.InitialBaserockPorosity", 0.01))
    , initBasalAquiferPorosity_(getParam<Scalar>("SpatialParams.InitialBasalAquiferPorosity", 0.1))
    , yBaserockMin_(getParam<Scalar>("SpatialParams.yBaserockMin", 370.0))
    , yBaserockMax_(getParam<Scalar>("SpatialParams.yBaserockMax", 470.0))
    , yReservoirMin_(getParam<Scalar>("SpatialParams.yReservoirMin", 470.0))
    , yReservoirMax_(getParam<Scalar>("SpatialParams.yReservoirMax", 500.0))
    , yCaprockMin_ (getParam<Scalar>("SpatialParams.yCaprockMin", 500.0))
    , yCaprockMax_ (getParam<Scalar>("SpatialParams.yCaprockMax", 600.0))
    {
        // Young's modulus [Pa]
       Scalar Er = 15e9;
       Scalar Ec = 5e9;
       Scalar Eu = 10e9;
       Scalar Eb = 20e9;
       // Poisson's ratio [-]
       Scalar nu = 0.3;
       // Lame parameters [Pa]
       ReservoirLameParams_.setLambda( (Er * nu) / ((1 + nu)*(1 - 2 * nu)) );
       ReservoirLameParams_.setMu( Er / (2 * (1 + nu)) );
       CaprockLameParams_.setLambda( (Ec * nu) / ((1 + nu)*(1 - 2 * nu)) );
       CaprockLameParams_.setMu( Ec / (2 * (1 + nu)) );
       BaserockLameParams_.setLambda( (Ec * nu) / ((1 + nu)*(1 - 2 * nu)) );
       BaserockLameParams_.setMu( Ec / (2 * (1 + nu)) );
       UpperAquiferLameParams_.setLambda( (Eu * nu) / ((1 + nu)*(1 - 2 * nu)) );
       UpperAquiferLameParams_.setMu( Eu / (2 * (1 + nu)) );
       BasalAquiferLameParams_.setLambda( (Eb * nu) / ((1 + nu)*(1 - 2 * nu)) );
       BasalAquiferLameParams_.setMu( Eb / (2 * (1 + nu)) );

    }

    //! Defines the Lame parameters.
    const LameParams& lameParamsAtPos(const GlobalPosition& globalPos) const
    {
        if (isUpperAquifer_ (globalPos))
        	return UpperAquiferLameParams_ ;
        else if (isCaprock_ (globalPos))
        	return CaprockLameParams_  ;
		else if (isReservoir_ (globalPos))
			return ReservoirLameParams_ ;
		else if (isBaserock_ (globalPos))
			return BaserockLameParams_ ;
		else
			return BasalAquiferLameParams_ ;

//    	return lameParams_;
    }

    //! Returns the porosity of the porous medium.
    template<class ElementSolution>
    Scalar porosity(const Element& element,
                    const SubControlVolume& scv,
                    const ElementSolution& elemSol) const

    {
        const GlobalPosition& globalPos = element.geometry().center();
        if (isUpperAquifer_ (globalPos))
        	return initUpperAquiferPorosity_ ;
        else if (isCaprock_ (globalPos))
        	return initCaprockPorosity_ ;
//        	return PorosityDeformation<Scalar>::evaluatePorosity(poroMechGridGeom, element, scv.center(), poroMechElemSol, initCaprockPorosity_);
		else if (isReservoir_ (globalPos))
////			return initReservoirPorosity_ ;
//			return PorosityDeformation<Scalar>::evaluatePorosity(poroMechGridGeom, element, scv.center(), poroMechElemSol, initReservoirPorosity_);
			return PorosityDeformation<Scalar>::evaluatePorosity(this->gridGeometry(), element, scv, elemSol, initReservoirPorosity_);
		else if (isBaserock_ (globalPos))
			return initBaserockPorosity_ ;
		else
			return initBasalAquiferPorosity_ ;
//			return PorosityDeformation<Scalar>::evaluatePorosity(poroMechGridGeom, element, scv.center(), poroMechElemSol, initBaserockPorosity_);

//        return PorosityDeformation<Scalar>::evaluatePorosity(this->gridGeometry(), element, scv, elemSol, initPorosity_);
    }

    //! Returns the Biot coefficient of the porous medium.
    Scalar biotCoefficientAtPos(const GlobalPosition& globalPos) const
    { return 1.0; }

private:
    Scalar initPorosity_;
    LameParams lameParams_;
    Scalar initCaprockPorosity_ , initReservoirPorosity_ , initBaserockPorosity_, initUpperAquiferPorosity_ , initBasalAquiferPorosity_;
    LameParams UpperAquiferLameParams_ , CaprockLameParams_, ReservoirLameParams_, BaserockLameParams_, BasalAquiferLameParams_;
    Scalar yReservoirMax_ , yReservoirMin_ , yCaprockMax_ , yCaprockMin_ , yBaserockMax_ , yBaserockMin_ ;
    Scalar yUpperAquifer_ = yCaprockMax_ ;
    Scalar yBasalAquifer_ = yBaserockMin_ ;

    bool isReservoir_ (const GlobalPosition &globalPos) const
    {
    	return globalPos[1] < yReservoirMax_ && globalPos[1] > yReservoirMin_ ;
    }

    bool isCaprock_(const GlobalPosition &globalPos) const
    {
    	return globalPos[1] < yCaprockMax_ && globalPos[1] > yCaprockMin_ ;
    }

    bool isBaserock_ (const GlobalPosition &globalPos) const
    {
    	return globalPos[1] < yBaserockMax_ && globalPos[1] > yBaserockMin_ ;
    }

    bool isUpperAquifer_ (const GlobalPosition& globalPos) const
    {
    	return globalPos[1] > yUpperAquifer_ ;
    }

    bool isBasalAquifer_ (const GlobalPosition& globalPos) const
    {
    	return globalPos[1] < yBasalAquifer_ ;
    }
};

} // end namespace Dumux

#endif
