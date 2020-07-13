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
 * \brief The spatial parameters class for the two-phase sub problem in the el2p test problem.
 */

#ifndef DUMUX_2P_TEST_SPATIALPARAMS_HH
#define DUMUX_2P_TEST_SPATIALPARAMS_HH

#include <dumux/discretization/elementsolution.hh>

#include <dumux/material/spatialparams/fv.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedvangenuchten.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>
#include <dumux/material/spatialparams/gstatrandomfield.hh>
#include <dumux/material/fluidmatrixinteractions/porositydeformation.hh>
#include <dumux/material/fluidmatrixinteractions/permeabilitykozenycarman.hh>

namespace Dumux {

/*!
 * \ingroup PoromechanicsTests
 * \brief The spatial parameters class for the two-phase sub problem in the el2p test problem.
 */
template<class GridGeometry, class Scalar, class CouplingManager>
class TwoPSpatialParams : public FVSpatialParams<GridGeometry, Scalar,
                                                 TwoPSpatialParams<GridGeometry, Scalar, CouplingManager>>
{
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using ThisType = TwoPSpatialParams<GridGeometry, Scalar, CouplingManager>;
    using ParentType = FVSpatialParams<GridGeometry, Scalar, ThisType>;

public:
//    using EffectiveLaw = RegularizedBrooksCorey<Scalar>;
    using EffectiveLaw = RegularizedVanGenuchten<Scalar>;
    using MaterialLaw = EffToAbsLaw<EffectiveLaw>;
    using MaterialLawParams = typename MaterialLaw::Params;
    // export permeability type
    using PermeabilityType = Scalar;

    TwoPSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry,
                      std::shared_ptr<CouplingManager> couplingManagerPtr)
    : ParentType(gridGeometry)
    , couplingManagerPtr_(couplingManagerPtr)
    , initPermeability_(getParam<Scalar>("SpatialParams.Permeability"))
    , initUpperAquiferPermeability_(getParam<Scalar>("SpatialParams.UpperAquiferPermeability", 1e-13))
    , initCaprockPermeability_(getParam<Scalar>("SpatialParams.CaprockPermeability", 1e-17))
    , initReservoirPermeability_(getParam<Scalar>("SpatialParams.ReservoirPermeability", 1.23e-13))
    , initBaserockPermeability_(getParam<Scalar>("SpatialParams.BaserockPermeability", 1e-18))
    , initBasalAquiferPermeability_(getParam<Scalar>("SpatialParams.BasalAquiferPermeability", 1e-14))
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
        // given Van Genuchten m
        Scalar m = 0.315;
        Scalar m_c = 0.457;
        Scalar m_u = 0.6;
        // Brooks Corey lambda
        using std::pow;
//        Scalar brooksCoreyLambda = m / (1 - m) * (1 - pow(0.5, 1/m));

        // residual saturations
        ReservoirMaterialParams_.setSwr(0.04);
        ReservoirMaterialParams_.setSnr(0.0);
        CaprockMaterialParams_.setSwr(0.3);
        CaprockMaterialParams_.setSnr(0.0);
        BaserockMaterialParams_.setSwr(0.3);
        BaserockMaterialParams_.setSnr(0.0);
        UpperAquiferMaterialParams_.setSwr(0.04);
        UpperAquiferMaterialParams_.setSnr(0.0);
        BasalAquiferMaterialParams_.setSwr(0.04);
        BasalAquiferMaterialParams_.setSnr(0.0);


        // parameters for the Brooks Corey law
//        myMaterialParams_.setPe(1.99e4);
//        myMaterialParams_.setLambda(brooksCoreyLambda);
        ReservoirMaterialParams_.setVgAlpha(1/0.0057e6);
		ReservoirMaterialParams_.setVgm(m);
		CaprockMaterialParams_.setVgAlpha(1/5e6);
		CaprockMaterialParams_.setVgm(m_c);
		BaserockMaterialParams_.setVgAlpha(1/5e6);
		BaserockMaterialParams_.setVgm(m_c);
		UpperAquiferMaterialParams_.setVgAlpha(1/0.01e6);
		UpperAquiferMaterialParams_.setVgm(m_u);
		BasalAquiferMaterialParams_.setVgAlpha(1/0.01e6);
		BasalAquiferMaterialParams_.setVgm(m_u);
    }

    //! Returns the porosity for a sub-control volume.
    template< class ElementSolution >
    Scalar porosity(const Element& element,
                    const SubControlVolume& scv,
                    const ElementSolution& elemSol) const
    {
        static constexpr auto poroMechId = CouplingManager::poroMechId;

        const auto& poroMechGridGeom = couplingManagerPtr_->problem(poroMechId).gridGeometry();
        const auto poroMechElemSol = elementSolution(element, couplingManagerPtr_->curSol()[poroMechId], poroMechGridGeom);

        // evaluate the deformation-dependent porosity at the scv center
        const GlobalPosition& globalPos = element.geometry().center();
        if (isUpperAquifer_ (globalPos))
        	return initUpperAquiferPorosity_ ;
        else if (isCaprock_ (globalPos))
        	return initCaprockPorosity_ ;
//        	return PorosityDeformation<Scalar>::evaluatePorosity(poroMechGridGeom, element, scv.center(), poroMechElemSol, initCaprockPorosity_);
		else if (isReservoir_ (globalPos))
////			return initReservoirPorosity_ ;
			return PorosityDeformation<Scalar>::evaluatePorosity(poroMechGridGeom, element, scv.center(), poroMechElemSol, initReservoirPorosity_);
		else if (isBaserock_ (globalPos))
			return initBaserockPorosity_ ;
		else
			return initBasalAquiferPorosity_ ;
//			return PorosityDeformation<Scalar>::evaluatePorosity(poroMechGridGeom, element, scv.center(), poroMechElemSol, initBaserockPorosity_);

//        return PorosityDeformation<Scalar>::evaluatePorosity(poroMechGridGeom, element, scv.center(), poroMechElemSol, initPorosity_);
    }

    //! Functions for defining the (intrinsic) permeability \f$[m^2]\f$.
     template< class ElementSolution >
     PermeabilityType permeability(const Element& element,
                                   const SubControlVolume& scv,
                                   const ElementSolution& elemSol) const
     {
         PermeabilityKozenyCarman<PermeabilityType> permLaw;
         const GlobalPosition& globalPos = element.geometry().center();
         if (isUpperAquifer_ (globalPos))
        	 return initUpperAquiferPermeability_ ;
         else if (isCaprock_ (globalPos))
        	 return initCaprockPermeability_ ;
////             return permLaw.evaluatePermeability(initCaprockPermeability_, initCaprockPorosity_, porosity(element, scv, elemSol));
 		else if (isReservoir_ (globalPos))
// 			return initReservoirPermeability_ ;
 	         return permLaw.evaluatePermeability(initReservoirPermeability_, initReservoirPorosity_, porosity(element, scv, elemSol));
 		else if (isBaserock_ (globalPos))
 			return initBaserockPermeability_ ;
 		else
 			return initBasalAquiferPermeability_ ;
//// 	         return permLaw.evaluatePermeability(initBaserockPermeability_, initBaserockPorosity_, porosity(element, scv, elemSol));

//         return permLaw.evaluatePermeability(initPermeability_, initPorosity_, porosity(element, scv, elemSol));
     }

    /*!
     * \brief Returns the parameter object for the Brooks-Corey material law.
     *
     * In this test, we use element-wise distributed material parameters.
     *
     * \param element The current element
     * \param scv The sub-control volume inside the element.
     * \param elemSol The solution at the dofs connected to the element.
     * \return The material parameters object
     */
    template<class ElementSolution>
    const MaterialLawParams& materialLawParams(const Element& element,
                                               const SubControlVolume& scv,
                                               const ElementSolution& elemSol) const
    {
        // do not use different parameters in the test with inverted wettability
//    	const GlobalPosition& globalPos = element.geometry().center();
//    	if (isUpperAquifer_ (globalPos))
//    		return UpperAquiferMaterialParams_ ;
//    	else if(isCaprock_ (globalPos))
//    		return CaprockMaterialParams_ ;
//    	else if (isReservoir_ (globalPos))
//    		return ReservoirMaterialParams_ ;
//    	else if (isBaserock_ (globalPos))
//    		return BaserockMaterialParams_ ;
//    	else
//    		return BasalAquiferMaterialParams_ ;
        return ReservoirMaterialParams_;
    }

    /*!
     * \brief Function for defining which phase is to be considered as the wetting phase.
     *
     * \return the wetting phase index
     * \param globalPos The global position
     */
    template<class FluidSystem>
    int wettingPhaseAtPos(const GlobalPosition& globalPos) const
    {
        return FluidSystem::phase0Idx;
    }

    //! Returns reference to the coupling manager.
    const CouplingManager& couplingManager() const
    { return *couplingManagerPtr_; }

private:
    std::shared_ptr<const CouplingManager> couplingManagerPtr_;
    Scalar initPermeability_;
    Scalar initPorosity_;
    Scalar yReservoirMax_ , yReservoirMin_ , yCaprockMax_ , yCaprockMin_ , yBaserockMax_ , yBaserockMin_ ;
    Scalar yUpperAquifer_ = yCaprockMax_ ;
    Scalar yBasalAquifer_ = yBaserockMin_ ;
    Scalar initCaprockPorosity_ , initReservoirPorosity_ , initBaserockPorosity_, initUpperAquiferPorosity_ , initBasalAquiferPorosity_;
    Scalar initCaprockPermeability_ , initReservoirPermeability_ , initBaserockPermeability_ , initUpperAquiferPermeability_ , initBasalAquiferPermeability_ ;

    MaterialLawParams myMaterialParams_;
    MaterialLawParams CaprockMaterialParams_ , ReservoirMaterialParams_ , BaserockMaterialParams_, UpperAquiferMaterialParams_ , BasalAquiferMaterialParams_;

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
