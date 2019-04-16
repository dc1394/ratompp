#include "stdafx.h"
#include "nonlinks.h"
#include <iostream>         // for std::cout
#include <boost/format.hpp> // for boost::format

namespace ks {
    //
    // Constructor
    //
    NonLinKs::NonLinKs(const char* path) : m_db(std::make_shared<ParamDb>(path))
    {
        m_db->ReadParams();
        m_db->WriteParams();
        
        m_pot = std::make_pair(std::make_shared<Pot<util::Spin::Alpha>>(m_db), std::make_shared<Pot<util::Spin::Beta>>(m_db));

        m_ss_alpha = std::make_shared<StateSet>(m_db);
        m_ss_beta = std::make_shared<StateSet>(m_db);
        
        m_ks = std::make_pair(std::make_shared< KohnSham<util::Spin::Alpha> >(m_db, m_ss_alpha), std::make_shared< KohnSham<util::Spin::Beta> >(m_db, std::move(m_ss_beta)));

        m_energy = std::make_unique<Energy>(m_pot, m_ss_alpha, m_db);

        m_rho = std::make_pair(std::make_shared<Rho>(m_db), std::make_shared<Rho>(m_db));
    }
    
    int NonLinKs::Run()
    {
        using namespace std::chrono;

        try
        {
            auto const beg = high_resolution_clock::now();
            Scf();
            auto const end = high_resolution_clock::now();

            std::cout << "\n\n\n********** WRITING RESULTS TO OUTPUT FILES **********\n";
            fflush(stdout);
            WriteRes(duration_cast<duration<double>>(end - beg));
        }
        catch (std::exception const & e)
        {
            std::cout << boost::format("\n\nERROR! %s\n\n\n") % e.what();
            return 1;
        }

        std::cout << "\n\n********** CALCULATIONS FINISHED SUCCESSFULLY! **********\n\n\n";
        return 0;
    }

    void NonLinKs::Scf(void)
    {
        auto const scfMaxIter = m_db->GetSize_t("Scf_MaxIter");
        auto pmix = std::make_pair(std::make_shared<RhoMix>(m_db), std::make_shared<RhoMix>(m_db));
        auto rhoOld = std::make_pair(std::make_shared<Rho>(m_db), std::make_shared<Rho>(m_db));
        auto iter = 1UL;
        
        m_rho.first->Init<util::Spin::Alpha>();
        m_rho.second->Init<util::Spin::Beta>();
        m_ks.first->Config(m_pot.first);
        m_ks.second->Config(m_pot.second);

        std::cout << "********************   S C F   L O O P   ********************\n";

        while (true)
        {
            std::cout << boost::format("*  SCF=%3lu   ") % static_cast<unsigned long>(iter);

            m_pot.first->SetRho(m_rho);
            m_pot.first->SolvePoisson();
            m_pot.second->SetRho(m_rho);
            m_pot.second->Hart() = m_pot.first->Hart;
            m_ks.first->Solve();
            m_ks.second->Solve();

            if (IsFinished())
                break;

            if (iter == scfMaxIter)
                break;
            {
                std::swap(m_rho, rhoOld);
            }

            //mix.SetRho(m_ks, rhoOld);
            pmix.first->SetRho(m_ks.first, rhoOld.first);
            pmix.second->SetRho(m_ks.second, rhoOld.second);
            //m_rho->Calc(&mix);
            m_rho.first->Calc(pmix.first);
            m_rho.second->Calc(pmix.second);

            iter++;
        }
        
        auto const node = m_rho.first->GetNode();
        m_energy->SetNode(node);
		
		m_pot.second->Write();
		m_pot.first->Write();

        std::cout << boost::format("*  SCF-ITERATIONS = %lu\n") % static_cast<unsigned long>(iter);
        std::cout << "***********   S C F   L O O P   F I N I S H E D   ***********\n";
    }



    //
    // 
    //
    bool NonLinKs::IsFinished() const
    {
        auto const scfEnerDiff = std::stof(m_db->Get("Scf_Diff"));
        static auto sumOld = 0.0; // This variable MUST BE "static"
        
        auto const sumNew = m_ss_alpha->EigenSum();
        auto const diff = std::fabs(sumNew - sumOld);
        sumOld = sumNew;

        std::cout << boost::format("EigenSum = %18.10lf    Diff = %18.10E\n") % sumNew % diff;
        fflush(stdout);

        return (diff < scfEnerDiff);
    }
        
    void NonLinKs::WriteRes(std::chrono::duration<double> const & sec) const
    {
        WriteInfo(sec);
        m_rho.first->Write();
        
        //m_ks->WriteEigen();
    }

    //
    // 
    //
    void NonLinKs::WriteInfo(std::chrono::duration<double> const & sec) const
    {
        auto out = m_db->OpenFile("out", "a");

        m_ss_alpha->WriteSates(stdout);
        m_energy->WriteEnergy(stdout);

        m_ss_alpha->WriteSates(out.get());
        m_energy->WriteEnergy(out.get());

        fprintf(out.get(), "\n\nC A L C U L A T I O N   T I M E :   %.3f  [s]\n", sec.count());
    }
}

