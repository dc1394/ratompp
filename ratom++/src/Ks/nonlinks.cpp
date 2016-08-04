#include "stdafx.h"
#include "nonlinks.h"

namespace ks {
    //
    // Constructor
    //
    NonLinKs::NonLinKs(const char* path) : m_db(path)
    {
        m_db.ReadParams();
        m_db.WriteParams();

        //m_pot = new Pot(&m_db);
        m_pot = std::make_pair(std::make_shared<Pot<util::Spin::Alpha>>(&m_db), std::make_shared<Pot<util::Spin::Alpha>>(&m_db));
        //assert(m_pot);

        m_ss = new StateSet(&m_db);
        assert(m_ss);

        // May 25th, 2014 Modified by dc1394
        //m_ks = new KohnSham(&m_db, m_ss);
        m_ks = std::make_pair(std::make_shared<KohnSham<util::Spin::Alpha>>(&m_db, m_ss), std::make_shared<KohnSham<util::Spin::Beta>>(&m_db, m_ss));
        //assert(m_ks);

        m_energy = new Energy(m_pot, m_ss, &m_db);
        assert(m_energy);

        // May 25th, 2014 Modified by dc1394
        //m_rho = new Rho(&m_db);
        m_rho = std::make_pair(std::make_shared<Rho>(&m_db), std::make_shared<Rho>(&m_db));
        //assert(m_rho);

    }

    //
    // Destructor
    //
    NonLinKs::~NonLinKs(void)
    {
        // May 25th, 2014 Modified by dc1394
        //delete m_pot;
        delete m_ss;
        //delete m_ks;
        delete m_energy;
        //delete m_rho;
    }

    //
    // Read in input parameters. Stars the calculation. Writes results to files.
    //
    int NonLinKs::Run()
    {
        time_t end, beg;

        try
        {
            time(&beg);
            Scf();
            time(&end);

            //printf("\n\n\n********** WRITING RESULTS TO OUTPUT FILES **********\n");
            //fflush(stdout);
            WriteRes(end - beg);
        }
        catch (std::exception& e)
        {
            printf("\n\nERROR! %s\n\n\n", e.what());
            return 1;
        }

        printf("\n\n********** CALCULATIONS FINISHED SUCCESSFULLY! **********\n\n\n");
        return 0;
    }

    //
    // Iterative solution of nonlinear Kohn-Sham equation, SCF loop
    //
    void NonLinKs::Scf(void)
    {
        const size_t scfMaxIter = m_db.GetSize_t("Scf_MaxIter");
        // May 23rd, 2014 Modified by dc1394
        //RhoMix mix(&m_db);
        std::pair<RhoMix, RhoMix> mix(&m_db, &m_db);
        //Rho *rhoOld;
        std::pair<std::shared_ptr<Rho>, std::shared_ptr<Rho>> rhoOld;
        size_t iter = 1;
        // May 23rd, 2014 Modified by dc1394
        rhoOld = std::make_pair(std::make_shared<Rho>(&m_db), std::make_shared<Rho>(&m_db));
        //assert(rhoOld);

        //m_rho->Init();
        m_rho.first->Init();
        m_rho.second->Init();
        m_ks.first->Config(m_pot.first);
        m_ks.second->Config(m_pot.second);

        printf("********************   S C F   L O O P   ********************\n");

        while (true)
        {
            printf("*  SCF=%3lu   ", static_cast<unsigned long>(iter));

            m_pot.first->SetRho(m_rho);
            m_pot.first->SolvePoisson();
            m_pot.second->SetRho(m_rho);
            m_pot.second->SolvePoisson();
            m_ks.first->Solve();
            m_ks.second->Solve();

            if (IsFinished())
                break;

            if (iter == scfMaxIter)
                break;
            {
                // May 23rd, 2014 Modified by dc1394
                //Rho *tmp = rhoOld;
                //rhoOld = m_rho;
                //m_rho = tmp;
                std::swap(m_rho, rhoOld);
            }

            //mix.SetRho(m_ks, rhoOld);
            mix.first.SetRho(m_ks.first, rhoOld.first);
            mix.second.SetRho(m_ks.second, rhoOld.second);
            //m_rho->Calc(&mix);
            m_rho.first->Calc(&mix.first);
            m_rho.second->Calc(&mix.second);

            iter++;
        }

        std::vector<double> node;
        m_rho->GetNode(node);
        m_energy->SetNode(node);
        // m_energy->WriteEnergy(stdout);
        // March 31st, 2014	Added by dc1394
        m_pot->Write();

        printf("*  SCF-ITERATIONS = %lu\n", static_cast<unsigned long>(iter));
        printf("***********   S C F   L O O P   F I N I S H E D   ***********\n");

        // May 23rd, 2014 Modified by dc1394
        //delete rhoOld;
    }



    //
    // Returns "true", if required accuarcy reached.
    //
    bool NonLinKs::IsFinished() const
    {
        const double scfEnerDiff = atof(m_db.Get("Scf_Diff"));
        static double sumOld = 0; // This variable MUST BE "static"
        double sumNew, diff;

        sumNew = m_ss->EigenSum();
        diff = ::fabs(sumNew - sumOld);
        sumOld = sumNew;

        printf("EigenSum = %18.10lf    Diff = %18.10E\n", sumNew, diff);
        fflush(stdout);

        return (diff < scfEnerDiff);
    }

    //
    // Write results into files
    //
    void NonLinKs::WriteRes(time_t sec) const
    {
        WriteInfo(sec);
        m_rho->Write();
        m_ks->WriteEigen();
    }

    //
    // Writes energy and eigenvalues into files.
    // "sec" - time of calcuylations in seconds.
    //
    void NonLinKs::WriteInfo(time_t sec) const
    {
        FILE* out;

        out = m_db.OpenFile("out", "a");

        m_ss->WriteSates(stdout);
        m_energy->WriteEnergy(stdout);

        m_ss->WriteSates(out);
        m_energy->WriteEnergy(out);

        fprintf(out, "\n\nC A L C U L A T I O N   T I M E :   %ld  [s]\n", sec);

        fclose(out);
    }
}
