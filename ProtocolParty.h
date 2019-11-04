#ifndef PROTOCOLPARTY_H_
#define PROTOCOLPARTY_H_

#include <stdlib.h>

#include <libscapi/include/primitives/Matrix.hpp>
#include <libscapi/include/cryptoInfra/Protocol.hpp>
#include <libscapi/include/circuits/ArithmeticCircuit.hpp>
#include <libscapi/include/infra/Measurement.hpp>
#include <vector>
#include <bitset>
#include <iostream>
#include <fstream>
#include <chrono>
#include <libscapi/include/primitives/Mersenne.hpp>
#include "ProtocolTimer.h"
#include <libscapi/include/comm/MPCCommunication.hpp>
#include <libscapi/include/infra/Common.hpp>
#include <libscapi/include/primitives/Prg.hpp>
#include "HashEncrypt.h"
#include <emmintrin.h>
#include <thread>

#define flag_print false
#define flag_print_timings true
#define flag_print_output true


using namespace std;
using namespace std::chrono;

template <class FieldType>
class ProtocolParty : public Protocol {

private:

    /**
     * N - number of parties
     * M - number of gates
     * T - number of malicious
     */

    bool isHonest = true;

    bool isPRF = true;

    int N, M, T, m_partyId;
    int times; //number of times to run the run function
    int iteration; //number of the current iteration

    Measurement* timer;
    VDM<FieldType> matrix_vand;
    TemplateField<FieldType> *field;
    vector<shared_ptr<ProtocolPartyData>>  parties;
    vector<FieldType> randomTAnd2TShares;
    vector<FieldType> secureRandomTShares;
    vector<FieldType> secureRandom2TShares;
    int randomOffset = 0;

    PrgFromOpenSSLAES prg;

    vector<PrgFromOpenSSLAES> my2TPrgs;
    vector<PrgFromOpenSSLAES> myTPrgs;
    vector<PrgFromOpenSSLAES> fromOthersTPrgs;
    vector<PrgFromOpenSSLAES> fromOthers2TPrgs;


    ProtocolTimer* protocolTimer;
    int currentCirciutLayer = 0;

    string s;
    int numOfInputGates, numOfOutputGates;
    string inputsFile, outputFile;
    vector<FieldType> beta;
    HIM<FieldType> matrix_for_interpolate;
    HIM<FieldType> matrix_for_interpolate_for_t_random_shares;
    HIM<FieldType> matrix_for_interpolate_for_2t_random_shares;
    HIM<FieldType> matrix_for_t;
    HIM<FieldType> matrix_for_2t;
    vector<FieldType> y_for_interpolate;

    vector<HIM<FieldType>> matrices_for_interpolate;

    HIM<FieldType> matrix_for_interpolate_mult2T;

    HIM<FieldType> matrix_him;

    VDMTranspose<FieldType> matrix_vand_transpose;

    HIM<FieldType> m;

    boost::asio::io_service io_service;
    ArithmeticCircuit circuit;
    vector<FieldType> gateValueArr; // the value of the gate (for my input and output gates)
    vector<FieldType> gateShareArr; // my share of the gate (for all gates)
    vector<FieldType> alpha; // N distinct non-zero field elements


    vector<FieldType> wSharesForVerification;
    vector<FieldType> multGatesSharesForVerification;
    vector<FieldType> randomElementsForVerification;

    vector<long> myInputs;

public:

    ProtocolParty(int argc, char* argv[]);


    void roundFunctionSync(vector<vector<byte>> &sendBufs, vector<vector<byte>> &recBufs, int round);
    void exchangeData(vector<vector<byte>> &sendBufs,vector<vector<byte>> &recBufs, int first, int last);

    void roundFunctionSyncElements(vector<vector<FieldType>> &sendBufs, vector<vector<FieldType>> &recBufs, int round);
    void exchangeDataElements(vector<vector<FieldType>> &sendBufs,vector<vector<FieldType>> &recBufs, int first, int last);


    int counter = 0;

    /**
     * This method runs the protocol:
     * 1. Preparation Phase
     * 2. Input Phase
     * 3. Computation Phase
     * 4. Verification Phase
     * 5. Output Phase
     */
    void run() override;

    bool hasOffline() {
        return true;
    }


    bool hasOnline() override {
        return true;
    }

    /**
     * This method runs the protocol:
     * Preparation Phase
     */
    void runOffline() override;

    /**
     * This method runs the protocol:
     * Input Phase
     * Computation Phase
     * Verification Phase
     * Output Phase
     */
    void runOnline() override;

    /**
     * This method reads text file and inits a vector of Inputs according to the file.
     */
    void readMyInputs();

    /**
     * We describe the protocol initialization.
     * In particular, some global variables are declared and initialized.
     */
    void initializationPhase();

    void initMatricesForRandomSharesAndKeys();


    bool preparationPhase();


    void offlineDNForMultiplication(int numOfTriples);


    /**
     * The input phase proceeds in two steps:
     * First, for each input gate, the party owning the input creates shares for that input by choosing a random coefficients for the polynomial
     * Then, all the shares are sent to the relevant party
     */
    void inputPhase();
    void generateRandom2TAndTShares(int numOfRandomPairs, vector<FieldType>& randomElementsToFill);
    void alternativeGenerateRandom2TAndTShares(int numOfRandomPairs, vector<FieldType>& randomElementsToFill);
    bool generateSecureDoubleSharings(int no_random);


    /**
     * Check whether given points lie on polynomial of degree d.
     * This check is performed by interpolating x on the first d + 1 positions of α and check the remaining positions.
     */
    bool checkConsistency(vector<FieldType>& x, int d);

    FieldType reconstructShare(vector<FieldType>& x, int d);

    void openShare(int numOfRandomShares, vector<FieldType> &Shares, vector<FieldType> &secrets, int degree);

    /**
     * Process all multiplications which are ready.
     * Return number of processed gates.
     */
    int processMultiplications(int lastMultGate);

    int processMultDN(int indexInRandomArray);

    int alternativeprocessMultDN(int indexInRandomArray);


    int processNotMult();

    /**
     * Walk through the circuit and evaluate the gates. Always take as many gates at once as possible,
     * i.e., all gates whose inputs are ready.
     * We first process all random gates, then alternately process addition and multiplication gates.
     */
    void computationPhase(HIM<FieldType> &m);

    /**
     * The cheap way: Create a HIM from the αi’s onto ZERO (this is actually a row vector), and multiply
     * this HIM with the given x-vector (this is actually a scalar product).
     * The first (and only) element of the output vector is the secret.
     */
    FieldType interpolate(vector<FieldType>& x);

    FieldType interpolateForParty(vector<FieldType>& x, int partyID);


    /**
     * Walk through the circuit and verify the multiplication gates.
     * We first generate the random elements using a common AES key that was generated by the parties,
     * run the relevant verification algorithm and return accept/reject according to the output
     * of the verification algorithm.
     */
    void verificationPhase();

    bool checkZero(int numOfChecks, vector<FieldType> &toCheck);



    vector<byte> generateCommonKey();
    void generatePseudoRandomElements(vector<byte> & aesKey, vector<FieldType> &randomElementsToFill, int numOfRandomElements);

    /**
     * Walk through the circuit and reconstruct output gates.
     */
    void outputPhase();

    ~ProtocolParty();

};


template <class FieldType>
ProtocolParty<FieldType>::ProtocolParty(int argc, char* argv[]) : Protocol("MPCHonestMajorityNoTriplesYehuda", argc, argv)
{
    string circuitFile = this->getParser().getValueByKey(arguments, "circuitFile");
    this->times = stoi(this->getParser().getValueByKey(arguments, "internalIterationsNumber"));
    string fieldType = this->getParser().getValueByKey(arguments, "fieldType");
    m_partyId = stoi(this->getParser().getValueByKey(arguments, "partyID"));
    int n = stoi(this->getParser().getValueByKey(arguments, "partiesNumber"));
    isHonest = stoi(this->getParser().getValueByKey(arguments, "isHonest"));
    isPRF = stoi(this->getParser().getValueByKey(arguments, "isPRF"));
    string outputTimerFileName = circuitFile + "Times" + to_string(m_partyId) + fieldType + ".csv";
    ProtocolTimer p(times, outputTimerFileName);

    this->protocolTimer = new ProtocolTimer(times, outputTimerFileName);

    if(isHonest==false){
        vector<string> subTaskNames{"Offline", "preparationPhase", "Online", "inputPhase", "ComputePhase",
                                    "VerificationPhase", "outputPhase"};
        timer = new Measurement(*this, subTaskNames);
    }
    else{
        vector<string> subTaskNames{"Offline", "preparationPhase", "Online", "inputPhase", "ComputePhase",
                                    "outputPhase"};
        timer = new Measurement(*this, subTaskNames);
    }

    if(fieldType.compare("ZpMersenne31") == 0) {
        field = new TemplateField<FieldType>(2147483647);
    } else if(fieldType.compare("ZpMersenne61") == 0) {
        field = new TemplateField<FieldType>(0);
    }else if(fieldType.compare("GF2_8LookupTable") == 0) {
        field = new TemplateField<FieldType>(0);
    }else if(fieldType.compare("ZpKaratsuba") == 0) {
        field = new TemplateField<FieldType>(0);
    } else if(fieldType.compare("GF2E") == 0) {
        field = new TemplateField<FieldType>(8);
    } else if(fieldType.compare("Zp") == 0) {
        field = new TemplateField<FieldType>(2147483647);
    }


    N = n;
    T = n/3 - 1;
    this->inputsFile = this->getParser().getValueByKey(arguments, "inputFile");
    this->outputFile = this->getParser().getValueByKey(arguments, "outputFile");
    if(n%3 > 0)
    {
        T++;
    }

    s = to_string(m_partyId);
    circuit.readCircuit(circuitFile.c_str());
    circuit.reArrangeCircuit();
    M = circuit.getNrOfGates();
    numOfInputGates = circuit.getNrOfInputGates();
    numOfOutputGates = circuit.getNrOfOutputGates();
    myInputs.resize(numOfInputGates);
    counter = 0;


    //comm->ConnectionToServer(s);

    //boost::asio::io_service io_service;

    MPCCommunication comm;
    string partiesFile = this->getParser().getValueByKey(arguments, "partiesFile");

    parties = comm.setCommunication(io_service, m_partyId, N, partiesFile);

    string tmp = "init times";
    //cout<<"before sending any data"<<endl;
    byte tmpBytes[20];
    for (int i=0; i<parties.size(); i++){
        if (parties[i]->getID() < m_partyId){
            parties[i]->getChannel()->write(tmp);
            parties[i]->getChannel()->read(tmpBytes, tmp.size());
        } else {
            parties[i]->getChannel()->read(tmpBytes, tmp.size());
            parties[i]->getChannel()->write(tmp);
        }
    }


    readMyInputs();

    auto t1 = high_resolution_clock::now();
    initializationPhase();

    auto t2 = high_resolution_clock::now();

    auto duration = duration_cast<milliseconds>(t2-t1).count();
    if(flag_print_timings) {
        cout << "time in milliseconds initializationPhase: " << duration << endl;
    }
}




template <class FieldType>
void ProtocolParty<FieldType>::readMyInputs()
{

    //cout<<"inputs file" << inputsFile<<endl;
    ifstream myfile;
    long input;
    int i =0;
    myfile.open(inputsFile);
    do {
        myfile >> input;
        myInputs[i] = input;
        i++;
    } while(!(myfile.eof()));
    myfile.close();


}

template <class FieldType>
void ProtocolParty<FieldType>::run() {

    for (iteration=0; iteration<times; iteration++){

        auto t1start = high_resolution_clock::now();
        timer->startSubTask("Offline", iteration);
        runOffline();
        timer->endSubTask("Offline", iteration);
        timer->startSubTask("Online", iteration);
        runOnline();
        timer->endSubTask("Online", iteration);

        auto t2end = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(t2end-t1start).count();
        protocolTimer->totalTimeArr[iteration] = duration;

        cout << "time in milliseconds for protocol: " << duration << endl;
    }


}

template <class FieldType>
void ProtocolParty<FieldType>::runOffline() {
    auto t1 = high_resolution_clock::now();
    timer->startSubTask("preparationPhase", iteration);
    if(preparationPhase() == false) {
        if(flag_print) {
            cout << "cheating!!!" << '\n';}
        return;
    }
    else {
        if(flag_print) {
            cout << "no cheating!!!" << '\n' << "finish Preparation Phase" << '\n';}
    }
    timer->endSubTask("preparationPhase", iteration);
    auto t2 = high_resolution_clock::now();

    auto duration = duration_cast<milliseconds>(t2-t1).count();
    if(flag_print_timings) {
        cout << "time in milliseconds preparationPhase: " << duration << endl;
    }
    protocolTimer->preparationPhaseArr[iteration] =duration;
}

template <class FieldType>
void ProtocolParty<FieldType>::runOnline() {

    auto t1 = high_resolution_clock::now();
    timer->startSubTask("inputPhase", iteration);
    inputPhase();
    timer->endSubTask("inputPhase", iteration);
    auto t2 = high_resolution_clock::now();

    auto duration = duration_cast<milliseconds>(t2-t1).count();
    protocolTimer->inputPreparationArr[iteration] = duration;
    if(flag_print_timings) {
        cout << "time in milliseconds inputPhase: " << duration << endl;
    }


    t1 = high_resolution_clock::now();
    timer->startSubTask("ComputePhase", iteration);
    computationPhase(m);
    timer->endSubTask("ComputePhase", iteration);
    t2 = high_resolution_clock::now();

    duration = duration_cast<milliseconds>(t2-t1).count();
    protocolTimer->computationPhaseArr[iteration] = duration;



    if(flag_print_timings) {
        cout << "time in milliseconds computationPhase: " << duration << endl;
    }

    if(isHonest==false) {
        t1 = high_resolution_clock::now();
        timer->startSubTask("VerificationPhase", iteration);
        verificationPhase();
        timer->endSubTask("VerificationPhase", iteration);
        t2 = high_resolution_clock::now();
        duration = duration_cast<milliseconds>(t2 - t1).count();
        protocolTimer->verificationPhaseArr[iteration] = duration;

        if (flag_print_timings) {
            cout << "time in milliseconds verificationPhase: " << duration << endl;
        }
    }

    t1 = high_resolution_clock::now();
    timer->startSubTask("outputPhase", iteration);
    outputPhase();
    timer->endSubTask("outputPhase", iteration);
    t2 = high_resolution_clock::now();

    duration = duration_cast<milliseconds>(t2-t1).count();
    protocolTimer->outputPhaseArr[iteration] = duration;

    if(flag_print_timings) {
        cout << "time in milliseconds outputPhase: " << duration << endl;
    }

}

template <class FieldType>
void ProtocolParty<FieldType>::computationPhase(HIM<FieldType> &m) {
    int count = 0;
    int countNumMult = 0;
    int countNumMultForThisLayer = 0;

    int numOfLayers = circuit.getLayers().size();
    for(int i=0; i<numOfLayers-1;i++){

        currentCirciutLayer = i;
        count = processNotMult();

        countNumMultForThisLayer = processMultiplications(countNumMult);//send the index of the current mult gate
        countNumMult += countNumMultForThisLayer;;
        count+=countNumMultForThisLayer;

    }
}

/**
 * the function implements the second step of Input Phase:
 * the party broadcasts for each input gate the difference between
 * the random secret and the actual input value.
 * @param diff
 */
template <class FieldType>
void ProtocolParty<FieldType>::inputPhase()
{
    int robin = 0;

    // the number of random double sharings we need altogether
    vector<FieldType> x1(N),y1(N);
    vector<vector<FieldType>> sendBufsElements(N);
    vector<vector<FieldType>> recBufElements(N);


    int index = 0;
    vector<int> sizes(N);

    // prepare the shares for the input
    for (int k = 0; k < numOfInputGates; k++)
    {
        if(circuit.getGates()[k].gateType == INPUT) {
            //get the expected sized from the other parties
            sizes[circuit.getGates()[k].party]++;

            if (circuit.getGates()[k].party == m_partyId) {
                auto input = myInputs[index];
                index++;
                if (flag_print) {
                    cout << "input  " << input << endl;
                }
                // the value of a_0 is the input of the party.
                x1[0] = field->GetElement(input);


                // generate random degree-T polynomial
                for(int i = 1; i < T+1; i++)
                {
                    // A random field element, uniform distribution
                    x1[i] = field->Random();

                }


                matrix_vand.MatrixMult(x1, y1, T+1); // eval poly at alpha-positions predefined to be alpha_i = i

                // prepare shares to be sent
                for(int i=0; i < N; i++)
                {
                    //cout << "y1[ " <<i<< "]" <<y1[i] << endl;
                    sendBufsElements[i].push_back(y1[i]);

                }
            }
        }
    }

    int fieldByteSize = field->getElementSizeInBytes();
    for(int i=0; i < N; i++)
    {
       recBufElements[i].resize(sizes[i]);

    }


    roundFunctionSyncElements(sendBufsElements, recBufElements,10);

    vector<int> counters(N);

    for(int i=0; i<N; i++){
        counters[i] =0;
    }

    for (int k = 0; k < numOfInputGates; k++)
    {
        if(circuit.getGates()[k].gateType == INPUT)
        {
            auto share = recBufElements[circuit.getGates()[k].party][counters[circuit.getGates()[k].party]];
            counters[circuit.getGates()[k].party] += 1;
            gateShareArr[circuit.getGates()[k].output] = share; // set the share sent from the party owning the input

        }
    }
    
}


template <class FieldType>
bool ProtocolParty<FieldType>::generateSecureDoubleSharings(int no_random)
{

    vector<vector<FieldType>> recBufsElements(N);
    vector<vector<FieldType>> sendBufsElements(N);

    int robin = 0;

    vector<FieldType> x1(N),x2(N),y1(N),y2(N);



    // the number of buckets (each bucket requires one double-sharing
    // from each party and gives N-2T random double-sharings)
    int no_buckets = (no_random / (N-2*T))+1;

    secureRandomTShares.resize(no_buckets*(N-2*T)); // my shares of the double-sharings
    secureRandom2TShares.resize(no_buckets*(N-2*T)); // my shares of the double-sharings

    for(int i=0; i < N; i++)
    {
        sendBufsElements[i].resize(no_buckets*2);
        recBufsElements[i].resize(no_buckets*2);
    }

    /**
     *  generate double sharings.
     *  first degree t.
     *  subsequent: degree 2t with same secret.
     */
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    for(int k=0; k < no_buckets; k++)
    {
        // generate random degree-T polynomial
        for(int i = 0; i < T+1; i++)
        {
            // A random field element, uniform distribution
            x1[i] = field->Random();

        }

        x2[0] = x1[0];


        for(int i = 1; i < 2*T+1; i++)
        {
            // otherwise random
            x2[i] = field->Random();
        }

        matrix_vand.MatrixMult(x1, y1, T+1); // eval poly at alpha-positions
        matrix_vand.MatrixMult(x2, y2, 2*T+1); // eval poly at alpha-positions

        // prepare shares to be sent
        for(int i=0; i < N; i++)
        {
            //cout << "y1[ " <<i<< "]" <<y1[i] << endl;
            sendBufsElements[i][2*k] = y1[i];
            sendBufsElements[i][2*k+1] = y2[i];
        }




    }//end print one



    if(flag_print) {
        for (int i = 0; i < N; i++) {
            for (int k = 0; k < sendBufsElements[0].size(); k++) {

                 cout << "before roundfunction4 send to " <<i <<" element: "<< k << " " << sendBufsElements[i][k] << endl;
            }
        }
    }

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>( t2 - t1 ).count();
    //cout << "generate random degree-T polynomial took : " <<duration<<" ms"<<endl;


    high_resolution_clock::time_point t3 = high_resolution_clock::now();

    roundFunctionSyncElements(sendBufsElements, recBufsElements,4);

    high_resolution_clock::time_point t4 = high_resolution_clock::now();
    auto duration2 = duration_cast<milliseconds>( t4 - t3 ).count();




    /**
     * Apply hyper-invertible matrix on each bucket.
     * From the resulting sharings, 2T are being reconstructed towards some party,
     * the remaining N-2T are kept as prepared sharings.
     * For balancing, we do round-robin the party how shall reconstruct and check!
     */


    for(int i=0; i<N; i++){
        sendBufsElements[i].clear();

    }

    int fieldBytesSize = field->getElementSizeInBytes();

    // x1 : used for the N degree-t sharings
    // x2 : used for the N degree-2t sharings
    for(int k=0; k < no_buckets; k++) {
        // generate random degree-T polynomial
        for (int i = 0; i < N; i++) {
            x1[i] = recBufsElements[i][2*k];
            x2[i] = recBufsElements[i][2*k+1];

        }
        matrix_him.MatrixMult(x1, y1);
        matrix_him.MatrixMult(x2, y2);
        // these shall be checked
        for (int i = 0; i < 2 * T; i++) {
            sendBufsElements[robin].push_back(y1[i]);
            sendBufsElements[robin].push_back(y2[i]);
            robin = (robin+1) % N; // next robin

        }
        // Y1 : the degree-t shares of my poly
        // Y2 : the degree 2t shares of my poly
        for (int i = 2 * T; i < N; i++) {

            secureRandomTShares[k*(N-2*T) + i - 2*T] = y1[i];
            secureRandom2TShares[k*(N-2*T) + i - 2*T] =  y2[i];
        }



        x2[0] = *(field->GetZero());
        x1[0] = *(field->GetZero());

    }

    for(int i=0; i < N; i++)
    {
        recBufsElements[i].resize(sendBufsElements[m_partyId].size());

    }


    t3 = high_resolution_clock::now();

    roundFunctionSyncElements(sendBufsElements, recBufsElements,5);

    t4 = high_resolution_clock::now();
    duration2 = duration_cast<milliseconds>( t4 - t3 ).count();


    if(flag_print) {
        cout << "after round" << endl;}
    int count = no_buckets * (2*T) / N; // nr of sharings *I* have to check
    // got one in the last round
    if(no_buckets * (2*T)%N > m_partyId) { // maybe -1
        count++;
    }


    for(int k=0; k < count; k++) {
        for (int i = 0; i < N; i++) {

            x1[i] = recBufsElements[i][2*k];
            x2[i] = recBufsElements[i][2*k+1];
        }


        vector<FieldType> x_until_d(N);
        for(int i=0; i<T; i++)
        {
            x_until_d[i] = x1[i];
        }
        for(int i=T; i<N; i++)
        {
            x_until_d[i] = *(field->GetZero());
        }
        if(flag_print) {

            //  cout <<"k "<<k<< "tinterpolate(x1).toString()  " << tinterpolate(x_until_d).toString() << endl;
            cout << "k " << k << "interpolate(x1).toString()  " << field->elementToString(interpolate(x1)) << endl;
            cout << "k " << k << "interpolate(x2).toString()  " << field->elementToString(interpolate(x2)) << endl;
        }
        // Check that x1 is t-consistent and x2 is 2t-consistent and secret is the same
        if(!checkConsistency(x1,T) || !checkConsistency(x2,2*T) ||
           (interpolate(x1)) != (interpolate(x2)))  {
            // cheating detected, abort
            if(flag_print) {
                cout << "k" << k<< endl;}
            return false;
        }
    }
    return true;
}
template <class FieldType>
void ProtocolParty<FieldType>::generateRandom2TAndTShares(int numOfRandomPairs, vector<FieldType>& randomElementsToFill){


    int index = 0;
    int robin = 0;
    int no_random = numOfRandomPairs;

    vector<FieldType> x1(N),y1(N), x2(N),y2(N), t1(N), r1(N), t2(N), r2(N);;

    vector<vector<FieldType>> sendBufsElements(N);
    vector<vector<FieldType>> recBufsElements(N);


    // the number of buckets (each bucket requires one double-sharing
    // from each party and gives N-2T random double-sharings)
    int no_buckets = (no_random / (N-T))+1;

    //maybe add some elements if a partial bucket is needed
    randomElementsToFill.resize(no_buckets*(N-T)*2);

    for(int i=0; i < N; i++)
    {
        sendBufsElements[i].resize(no_buckets*2);
        recBufsElements[i].resize(no_buckets*2);
    }

    /**
     *  generate random sharings.
     *  first degree t.
     *
     */
    for(int k=0; k < no_buckets; k++)
    {
        // generate random degree-T polynomial
        for(int i = 0; i < T+1; i++)
        {
            // A random field element, uniform distribution, note that x1[0] is the secret which is also random
            x1[i] = field->Random();

        }

        matrix_vand.MatrixMult(x1, y1,T+1); // eval poly at alpha-positions

        x2[0] = x1[0];
        // generate random degree-T polynomial
        for(int i = 1; i < 2*T+1; i++)
        {
            // A random field element, uniform distribution, note that x1[0] is the secret which is also random
            x2[i] = field->Random();

        }

        matrix_vand.MatrixMult(x2, y2,2*T+1);

        // prepare shares to be sent
        for(int i=0; i < N; i++)
        {
            //cout << "y1[ " <<i<< "]" <<y1[i] << endl;
            sendBufsElements[i][2*k] = y1[i];
            sendBufsElements[i][2*k + 1] = y2[i];

        }
    }

    if(flag_print) {
        for (int i = 0; i < N; i++) {
            for (int k = 0; k < sendBufsElements[0].size(); k++) {

                // cout << "before roundfunction4 send to " <<i <<" element: "<< k << " " << sendBufsElements[i][k] << endl;
            }
        }
        cout << "sendBufs" << endl;
        cout << "N" << N << endl;
        cout << "T" << T << endl;
    }

    roundFunctionSyncElements(sendBufsElements, recBufsElements,4);


    for(int k=0; k < no_buckets; k++) {
        for (int i = 0; i < N; i++) {
            t1[i] = recBufsElements[i][2*k];//field->bytesToElement(recBufsBytes[i].data() + (2*k * fieldByteSize));
            t2[i] = recBufsElements[i][2*k+1];//field->bytesToElement(recBufsBytes[i].data() + ((2*k +1) * fieldByteSize));

        }
        matrix_vand_transpose.MatrixMult(t1, r1,N-T);
        matrix_vand_transpose.MatrixMult(t2, r2,N-T);

        //copy the resulting vector to the array of randoms
        for (int i = 0; i < (N - T); i++) {

            randomElementsToFill[index*2] = r1[i];
            randomElementsToFill[index*2 +1] = r2[i];
            index++;

        }
    }

}


template <class FieldType>
void ProtocolParty<FieldType>::alternativeGenerateRandom2TAndTShares(int numOfRandomPairs, vector<FieldType>& randomElementsToFill) {


    int countGetBytes = 0;
    int index = 0;
    vector<vector<FieldType>>recBufsElements(N);

    int robin = 0;
    int no_random = numOfRandomPairs;

    vector<FieldType> x1(N), y1(N), x2(N), t1(N), r1(N), t2(N), r2(N);

    vector<vector<FieldType>> sendBufsElements(N);

    vector<byte> send2TShares;
    vector<byte> rec2TShares;



    // the number of buckets (each bucket requires one double-sharing
    // from each party and gives N-2T random double-sharings)
    int no_buckets = (no_random / (N - T)) + 1;

    //cout<<"num of buckets is : "<<no_buckets<<endl;

    vector<vector<FieldType>> y2(N);

     //maybe add some elements if a partial bucket is needed
    randomElementsToFill.resize(no_buckets * (N - T) * 2);



    //generate the shares for the t parties that have keys

    int counter = 1;
    int counter2T = 1;

    vector<FieldType> forInterpolate(T + 1);
    vector<FieldType> forInterpolate2T(2 * T + 1);
    vector<FieldType> calcShares(N - T);
    vector<FieldType> calcShares2T(1);


    for (int i = 0; i < N; i++) {
        //resize sendbufs according to the T parties my party sends data to
        if (!((m_partyId - i <= T && m_partyId - i > 0) || ((N + m_partyId) - i) <= T)) {
            sendBufsElements[i].resize(no_buckets);
        } else {
            sendBufsElements[i].resize(0);
        }


        //resize recdbufs according to the parties that send data to me
        if (!((i - m_partyId <= T && i - m_partyId > 0) || ((N + i) - m_partyId) <= T)) {
            recBufsElements[i].resize(no_buckets);
        } else {
            recBufsElements[i].resize(0);
        }

    }

    for (int i = 0; i < N; i++) {
        y2[i].resize(no_buckets);
    }


    for (int k = 0; k < no_buckets; k++) {

        //choose the share for 0 - the secret
        forInterpolate2T[0] = forInterpolate[0] = field->bytesToElement(
                myTPrgs[m_partyId].getPRGBytesEX(field->getElementSizeInBytes()));
        for (int i = 0; i < N; i++) {

            //should send seed for T also
            if ((m_partyId - i <= T && m_partyId - i > 0) || ((N + m_partyId) - i) <= T) {

                //set the shares from the prg as set by the t parties that have the prg that I had sent them
                forInterpolate[counter] = y1[i] =
                        field->bytesToElement(myTPrgs[i].getPRGBytesEX(field->getElementSizeInBytes()));

//                if (k == 10)
//                    cout << m_partyId << " my prg outputs is: " << forInterpolate[counter] << " for party " << i
//                         << endl;


                countGetBytes++;
                counter++;
            }

            if (m_partyId > 2 * T) {
                if (i <= 2 * T && i != m_partyId - (2 * T + 1)) {//do not sent to this party, calc for it the share
                    //we do not want to send always to party 2T, incase this party id is less than 2T.

                    //set the shares from the prg as set by the t parties that have the prg that I had sent them
                    forInterpolate2T[counter2T] = y2[i][k] = field->bytesToElement(
                            my2TPrgs[i].getPRGBytesEX(field->getElementSizeInBytes()));

                    counter2T++;

                }

            } else {//m_partyId<<2*T

                //calc my own share using interpolate. do not set my own share

                if (i <= 2 * T && i != m_partyId) {


                    //set the shares from the prg as set by the t parties that have the prg that I had sent them
                    forInterpolate2T[counter2T] = y2[i][k] = field->bytesToElement(
                            my2TPrgs[i].getPRGBytesEX(field->getElementSizeInBytes()));

                    counter2T++;

                }
            }
        }

        //generate the shares for all other parties
        matrix_for_interpolate_for_t_random_shares.MatrixMult(forInterpolate, calcShares);


        counter = 0;
        for (int i = 0; i < N; i++) {

            //fill the calculated shares for the parties that do not have keys including me
            if (!((m_partyId - i <= T && m_partyId - i > 0) || ((N + m_partyId) - i) <= T)) {

                y1[i] = calcShares[counter];

                sendBufsElements[i][k] = y1[i];
                counter++;
            }


        }
        counter = 1;
        counter2T = 1;

        matrix_for_interpolate_for_2t_random_shares.MatrixMult(forInterpolate2T, calcShares2T);


        if (m_partyId <= 2 * T) {
            //set my share to the calculated one
            y2[m_partyId][k] = calcShares2T[0];
        } else {

            //will need to send this to the relevant party
            y2[m_partyId - (2 * T + 1)][k] = calcShares2T[0];
        }
    }

    int fieldByteSize = field->getElementSizeInBytes();

    //takes care of t-degree shares
    roundFunctionSyncElements(sendBufsElements, recBufsElements, 4);

    //now handle the 2T shares

    if (m_partyId > 2 * T) {
        //only need to send

        send2TShares.resize(no_buckets * field->getElementSizeInBytes());

        field->elementVectorToByteVector(y2[m_partyId - (2 * T + 1)], send2TShares);

        parties[m_partyId - (2 * T + 1)]->getChannel()->write(send2TShares.data(), send2TShares.size());

    } else {

        if (m_partyId + (2 * T + 1) < N) {
            rec2TShares.resize(no_buckets * field->getElementSizeInBytes());

            //receive shares from the other party. Need to reduce one since my id is less and does not exist in the parties arrays
            parties[m_partyId + (2 * T + 1) - 1]->getChannel()->read(rec2TShares.data(), rec2TShares.size());
            //cout<<"read the data:: my Id = " << m_partyId-1<< "other ID = "<< parties[i]->getID()<<endl;
        }
        //only receive
    }



    for (int k = 0; k < no_buckets; k++) {
        for (int i = 0; i < N; i++) {

            if ((i - m_partyId <= T && i - m_partyId > 0) || ((N + i) - m_partyId) <= T) {
                t1[i] = field->bytesToElement(fromOthersTPrgs[i].getPRGBytesEX(field->getElementSizeInBytes()));

            } else {

                t1[i] = recBufsElements[i][k];
            }


            if (m_partyId <= 2 * T) {//calc only for these parties
                //set the 2T shares

                if ((m_partyId == (i - (2 * T + 1)))) {
                    t2[i] = field->bytesToElement(rec2TShares.data() + (k * fieldByteSize));
                } else if (i == m_partyId) {
                    t2[i] = y2[m_partyId][k];

                } else {

                    t2[i] = field->bytesToElement(fromOthers2TPrgs[i].getPRGBytesEX(field->getElementSizeInBytes()));

                }

            }
        }

        matrix_vand_transpose.MatrixMult(t1, r1, N - T);


        if (m_partyId <= 2 * T) {
            //do the actual random elements calculations only for parties that need the 2T random share
            //for the multiplication
            matrix_vand_transpose.MatrixMult(t2, r2, N - T);
        }

        //copy the resulting vector to the array of randoms
        for (int i = 0; i < (N - T); i++) {



            randomElementsToFill[index * 2] = r1[i];
            randomElementsToFill[index * 2 + 1] = r2[i];
            index++;


        }
    }
}





/**
 * some global variables are initialized
 * @param GateValueArr
 * @param GateShareArr
 * @param matrix_him
 * @param matrix_vand
 * @param alpha
 */
template <class FieldType>
void ProtocolParty<FieldType>::initializationPhase()
{
    beta.resize(1);
    y_for_interpolate.resize(N);
    gateShareArr.resize((M - circuit.getNrOfOutputGates())); // my share of the gate (for all gates)
    alpha.resize(N); // N distinct non-zero field elements
    vector<FieldType> alpha1(N-T);
    vector<FieldType> alpha2(T);

    beta[0] = field->GetElement(0); // zero of the field

    matrix_for_interpolate.allocate(1,N, field);


    my2TPrgs.resize(N);
    myTPrgs.resize(N);
    fromOthersTPrgs.resize(N);
    fromOthers2TPrgs.resize(N);


    matrix_him.allocate(N,N,field);
    matrix_vand.allocate(N,N,field);
    matrix_vand_transpose.allocate(N,N,field);
    m.allocate(T, N-T,field);

    // Compute Vandermonde matrix VDM[i,k] = alpha[i]^k
    matrix_vand.InitVDM();
    matrix_vand_transpose.InitVDMTranspose();

    // Prepare an N-by-N hyper-invertible matrix
    matrix_him.InitHIM();

    // N distinct non-zero field elements
    for(int i=0; i<N; i++)
    {
        alpha[i]=field->GetElement(i+1);
    }

    for(int i = 0; i < N-T; i++)
    {
        alpha1[i] = alpha[i];
    }
    for(int i = N-T; i < N; i++)
    {
        alpha2[i - (N-T)] = alpha[i];
    }

    m.InitHIMByVectors(alpha1, alpha2);



    vector<FieldType> alpha_until_t(T + 1);
    vector<FieldType> alpha_from_t(N - 1 - T);

    // Interpolate first d+1 positions of (alpha,x)
    matrix_for_t.allocate(N - 1 - T, T + 1, field); // slices, only positions from 0..d
    //matrix_for_t.InitHIMByVectors(alpha_until_t, alpha_from_t);
    matrix_for_t.InitHIMVectorAndsizes(alpha, T+1, N-T-1);

    vector<FieldType> alpha_until_2t(2*T + 1);
    vector<FieldType> alpha_from_2t(N - 1 - 2*T);

    // Interpolate first d+1 positions of (alpha,x)
    matrix_for_2t.allocate(N - 1 - 2*T, 2*T + 1, field); // slices, only positions from 0..d
    //matrix_for_2t.InitHIMByVectors(alpha_until_2t, alpha_from_2t);
    matrix_for_2t.InitHIMVectorAndsizes(alpha, 2*T + 1, N-(2*T +1));


    if(flag_print){
        cout<< "matrix_for_t : " <<endl;
        matrix_for_t.Print();

        cout<< "matrix_for_2t : " <<endl;
        matrix_for_2t.Print();

    }

    int delta =   (5 + field->getElementSizeInBytes() - 1) / field->getElementSizeInBytes();
    randomElementsForVerification.resize(circuit.getNrOfMultiplicationGates()*delta +
                                         (circuit.getNrOfMultiplicationGates() + circuit.getNrOfInputGates())*delta);

    wSharesForVerification.resize(circuit.getNrOfInputGates());
    multGatesSharesForVerification.resize(circuit.getNrOfMultiplicationGates()*2);



    matrix_for_interpolate.InitHIMByVectors(alpha, beta);



    vector<FieldType> currAlpha = alpha;
    currAlpha.resize(2*T + 1);
    matrix_for_interpolate_mult2T.allocate(1,2*T+1, field);
    matrix_for_interpolate_mult2T.InitHIMByVectors(currAlpha, beta);

    initMatricesForRandomSharesAndKeys();

}

template <class FieldType>
void ProtocolParty<FieldType>::initMatricesForRandomSharesAndKeys(){



    auto randomKey = prg.generateKey(128);
    prg.setKey(randomKey);
    int sizeOfSeed = 16;
    vector<vector<byte>> sendBufsBytes(N);
    vector<vector<byte>> recBufsBytes(N);

    vector<FieldType> currAlphaForT(T+1);
    vector<FieldType> currAlphaFor2T(2*T+1);

    vector<FieldType> currBetaForT(N-T);
    vector<FieldType> currBetaFor2T(1);


    int counterAlphaForT = 1;
    int counterAlphaFor2T = 1;
    int counterBetaForT = 0;
    int counterBetaFor2T = 0;

    currAlphaForT[0] = *field->GetZero();
    currAlphaFor2T[0] = *field->GetZero();


    //handle T prg's

    vector<byte> mykey(sizeOfSeed);
    prg.getPRGBytes(mykey, 0, sizeOfSeed);

    SecretKey key(mykey, "aes");
    my2TPrgs[m_partyId].setKey(key);

    prg.getPRGBytes(mykey, 0, sizeOfSeed);

    SecretKey key2(mykey, "aes");

    myTPrgs[m_partyId].setKey(key2);


    for(int i=0; i<N; i++) {

        //should send seed for T also
        if((m_partyId-i<=T && m_partyId-i>0) || ((N+m_partyId) - i)<=T){

            //fill the send buf with seed to 2T other parties and me
            sendBufsBytes[i].resize(sizeOfSeed);
            prg.getPRGBytes(sendBufsBytes[i], 0, sizeOfSeed);

            SecretKey key(sendBufsBytes[i], "aes");
            myTPrgs[i].setKey(key);


            //set the alpha values according to the parties that hold the t-degree key
            currAlphaForT[counterAlphaForT] = alpha[i];

            counterAlphaForT++;
        }
        else{

            //set relevant alpha

            //create beta for every point we will need to interpolate, this is for every party that does not
            //generate the point using the t-degree key
            currBetaForT[counterBetaForT] = alpha[i];

            counterBetaForT++;


        }
    }

    for(int i=0;i<N;i++){

        if((i-m_partyId<=T && i-m_partyId>0) || ((N+i) - m_partyId)<=T){
            recBufsBytes[i].resize(sizeOfSeed);
        }

        else{
            recBufsBytes[i].resize(0);
        }

    }

    //exchage T keys between all parties
    roundFunctionSync(sendBufsBytes , recBufsBytes,23);


    //set the prg that I get from other parties
    for(int i=0; i<N; i++){

        if((i-m_partyId<=T && i-m_partyId>0) || ((N+i) - m_partyId)<=T){

            SecretKey key(recBufsBytes[i], "aes");
            fromOthersTPrgs[i].setKey(key);

        }
    }



    for(int i=0; i<N; i++) {


        //sent seed for 2T shares
        if(m_partyId<=2*T) {

            if (i <= 2*T && i!=m_partyId) {


                //fill the send buf with seed to 2T other parties and me
                sendBufsBytes[i].resize(sizeOfSeed);
                prg.getPRGBytes(sendBufsBytes[i], 0, sizeOfSeed);

                SecretKey key(sendBufsBytes[i], "aes");
                my2TPrgs[i].setKey(key);


                currAlphaFor2T[counterAlphaFor2T] = alpha[i];

                counterAlphaFor2T++;
            } else {
                sendBufsBytes[i].resize(0);

                //currBetaFor2T[counterBetaFor2T] = alpha[i];

                //counterBetaFor2T++;
            }
        }

        if(m_partyId>2*T) {
            if (i <= 2*T ) {//do not sent to this party, calc for it the share
                                                     //we do not want to send always to party 2T, incase this party id is less than 2T.

                //fill the send buf with seed to 2T other parties and me
                sendBufsBytes[i].resize(sizeOfSeed);
                prg.getPRGBytes(sendBufsBytes[i], 0, sizeOfSeed);

                SecretKey key(sendBufsBytes[i], "aes");
                my2TPrgs[i].setKey(key);



                if( i==(m_partyId - (2*T+1))) {

                    sendBufsBytes[i].resize(0);

                }
                else{
                    currAlphaFor2T[counterAlphaFor2T] = alpha[i];

                    counterAlphaFor2T++;
                }



            } else {
                sendBufsBytes[i].resize(0);
            }
        }
    }

    //if I am one of the first 2T parties I will recieve keys, otherwise I will get no keys
    if(m_partyId<=2*T) {
        currBetaFor2T[0] = alpha[m_partyId];
        for (int i = 0; i < N; i++) {

            recBufsBytes[i].resize(sizeOfSeed);

            if(m_partyId==(i - (2*T+1)))//do not recieve from this party, the party will calculate the share for me and that share wiill be recieved
                recBufsBytes[i].resize(0);
        }
    }
    else{
        currBetaFor2T[0] = alpha[m_partyId - (2*T+1)];
        for (int i = 0; i < N; i++) {

            recBufsBytes[i].resize(0);
        }

    }

    //exchage 2T keys between all parties
    roundFunctionSync(sendBufsBytes , recBufsBytes,23);


    //set the prg that I get from other parties
    if(m_partyId<=2*T) {
        for (int i = 0; i < N; i++) {


            if(m_partyId!=(i - (2*T+1))&& i!=m_partyId) {//do not recieve from this party
                SecretKey key(recBufsBytes[i], "aes");
                fromOthers2TPrgs[i].setKey(key);
            }
        }
    }

    //finally create the matrix for the intrpolation of N-T-1 points. Every row in the
    //matrix will generate a N-T-1 shares for the parties that do not have the prg
    //to generate a random share by itself
    matrix_for_interpolate_for_t_random_shares.allocate(N-T,T+1, field);
    matrix_for_interpolate_for_t_random_shares.InitHIMByVectors(currAlphaForT, currBetaForT);


    matrix_for_interpolate_for_2t_random_shares.allocate(1, 2*T+1, field);
    matrix_for_interpolate_for_2t_random_shares.InitHIMByVectors(currAlphaFor2T, currBetaFor2T);




}

template <class FieldType>
bool ProtocolParty<FieldType>::preparationPhase()
{
    int delta =   (5 + field->getElementSizeInBytes() - 1) / field->getElementSizeInBytes();
    int keysize = 16/field->getElementSizeInBytes() + 1;

    randomOffset = 0;

    if(isHonest==false)
        generateSecureDoubleSharings(3*delta + keysize);

    //run offline for all the future multiplications including the multiplication of the protocol
    offlineDNForMultiplication(circuit.getNrOfMultiplicationGates());

    return true;
}


/**
 * Check whether given points lie on polynomial of degree d. This check is performed by interpolating x on
 * the first d + 1 positions of α and check the remaining positions.
 */
template <class FieldType>
bool ProtocolParty<FieldType>::checkConsistency(vector<FieldType>& x, int d)
{
    if(d == T)
    {
        vector<FieldType> y(N - 1 - d); // the result of multiplication
        vector<FieldType> x_until_t(T + 1);

        for (int i = 0; i < T + 1; i++) {
            x_until_t[i] = x[i];
        }


        matrix_for_t.MatrixMult(x_until_t, y);

        // compare that the result is equal to the according positions in x
        for (int i = 0; i < N - d - 1; i++)   // n-d-2 or n-d-1 ??
        {
            if ((y[i]) != (x[d + 1 + i])) {
                return false;
            }
        }
        return true;
    } else if (d == 2*T)
    {
        vector<FieldType> y(N - 1 - d); // the result of multiplication

        vector<FieldType> x_until_2t(2*T + 1);

        for (int i = 0; i < 2*T + 1; i++) {
            x_until_2t[i] = x[i];
        }

        matrix_for_2t.MatrixMult(x_until_2t, y);

        // compare that the result is equal to the according positions in x
        for (int i = 0; i < N - d - 1; i++)   // n-d-2 or n-d-1 ??
        {
            if ((y[i]) != (x[d + 1 + i])) {
                return false;
            }
        }
        return true;

    } else {
        vector<FieldType> alpha_until_d(d + 1);
        vector<FieldType> alpha_from_d(N - 1 - d);
        vector<FieldType> x_until_d(d + 1);
        vector<FieldType> y(N - 1 - d); // the result of multiplication

        for (int i = 0; i < d + 1; i++) {
            alpha_until_d[i] = alpha[i];
            x_until_d[i] = x[i];
        }
        for (int i = d + 1; i < N; i++) {
            alpha_from_d[i - (d + 1)] = alpha[i];
        }
        // Interpolate first d+1 positions of (alpha,x)
        HIM<FieldType> matrix(N - 1 - d, d + 1, field); // slices, only positions from 0..d
        matrix.InitHIMByVectors(alpha_until_d, alpha_from_d);
        matrix.MatrixMult(x_until_d, y);

        // compare that the result is equal to the according positions in x
        for (int i = 0; i < N - d - 1; i++)   // n-d-2 or n-d-1 ??
        {
            if (y[i] != x[d + 1 + i]) {
                return false;
            }
        }
        return true;
    }
    return true;
}

// Interpolate polynomial at position Zero
template <class FieldType>
FieldType ProtocolParty<FieldType>::interpolate(vector<FieldType>& x)
{
    //vector<FieldType> y(N); // result of interpolate
    matrix_for_interpolate.MatrixMult(x, y_for_interpolate);
    return y_for_interpolate[0];
}


template <class FieldType>
FieldType ProtocolParty<FieldType>::interpolateForParty(vector<FieldType>& x, int partyID)
{
    //vector<FieldType> y(N); // result of interpolate
    matrix_for_interpolate_mult2T.MatrixMult(x, y_for_interpolate);
    return y_for_interpolate[0];
}


template <class FieldType>
FieldType ProtocolParty<FieldType>::reconstructShare(vector<FieldType>& x, int d){

    if (!checkConsistency(x, d))
    {
        // someone cheated!

            cout << "cheating!!!" << '\n';
        exit(0);
    }
    else
        return interpolate(x);
}


template <class FieldType>
int ProtocolParty<FieldType>::processNotMult(){
    int count=0;
    for(int k=circuit.getLayers()[currentCirciutLayer]; k < circuit.getLayers()[currentCirciutLayer+1]; k++)
    {


        // add gate
        if(circuit.getGates()[k].gateType == ADD)
        {
            gateShareArr[circuit.getGates()[k].output] = gateShareArr[circuit.getGates()[k].input1] + gateShareArr[circuit.getGates()[k].input2];

            count++;
        }

        else if(circuit.getGates()[k].gateType == SUB)//sub gate
        {
            gateShareArr[circuit.getGates()[k].output] = gateShareArr[circuit.getGates()[k].input1] - gateShareArr[circuit.getGates()[k].input2];

            count++;
        }
        else if(circuit.getGates()[k].gateType == SCALAR)
        {
            long scalar(circuit.getGates()[k].input2);
            FieldType e = field->GetElement(scalar);
            gateShareArr[circuit.getGates()[k].output] = gateShareArr[circuit.getGates()[k].input1] * e;

            count++;
        }
        else if(circuit.getGates()[k].gateType == SCALAR_ADD)
        {
            long scalar(circuit.getGates()[k].input2);
            FieldType e = field->GetElement(scalar);
            gateShareArr[circuit.getGates()[k].output] = gateShareArr[circuit.getGates()[k].input1] + e;



            count++;
        }

         

    }
    return count;

}

/**
 * the Function process all multiplications which are ready.
 * @return the number of processed gates.
 */
template <class FieldType>
int ProtocolParty<FieldType>::processMultiplications(int lastMultGate)
{

    if(isPRF)
        return alternativeprocessMultDN(lastMultGate);
    else
        return processMultDN(lastMultGate);

}


template <class FieldType>
int ProtocolParty<FieldType>::alternativeprocessMultDN(int indexInRandomArray) {

    int index = 0;
    int fieldByteSize = field->getElementSizeInBytes();
    int maxNumberOfLayerMult = circuit.getLayers()[currentCirciutLayer + 1] - circuit.getLayers()[currentCirciutLayer];
    vector<FieldType> xyMinusRShares(maxNumberOfLayerMult*2);//hold both in the same vector to send in one batch

    vector<FieldType> xyMinusR;//hold both in the same vector to send in one batch
    vector<byte> xyMinusRBytes;


    vector<vector<FieldType>> sendBufsElements(N);
    vector<vector<FieldType>> recBufsElements(N);

    //generate the shares for x+a and y+b. do it in the same array to send once
    for (int k = circuit.getLayers()[currentCirciutLayer];
         k < circuit.getLayers()[currentCirciutLayer + 1]; k++)//go over only the logit gates
    {
        auto gate = circuit.getGates()[k];

        if (gate.gateType == MULT) {

            if(m_partyId<=2*T) {//only the first 2t+1 paties have the 2T share


                //compute the share of xy-r
                xyMinusRShares[index] = gateShareArr[gate.input1] * gateShareArr[gate.input2] -
                                        randomTAnd2TShares[2 * indexInRandomArray + 1];
            }

            indexInRandomArray++;


            index++;
        }
    }

    if(index==0)
        return 0;

    //set the acctual number of mult gate proccessed in this layer
    int acctualNumOfMultGates = index;
    int numOfElementsForParties = acctualNumOfMultGates/N;
    int indexForDecreasingSize = acctualNumOfMultGates - numOfElementsForParties *N;

    int counter=0;
    int currentNumOfElements;
    for(int i=0; i<N; i++){

        currentNumOfElements = numOfElementsForParties;
        if(i<indexForDecreasingSize)
            currentNumOfElements++;

        //should send
        if(m_partyId<=2*T ){


            //fill the send buf according to the number of elements to send to each party
            sendBufsElements[i].resize(currentNumOfElements);

            for(int j=0; j<currentNumOfElements; j++) {

                sendBufsElements[i][j] = xyMinusRShares[counter];
                counter++;

            }
            //field->elementVectorToByteVector(sendBufsElements[i], sendBufsBytes[i]);
        }
        else{
            sendBufsElements[i].resize(0);
        }
    }

    //resize the recbuf array.
    int myNumOfElementsToExpect = numOfElementsForParties;
    if (m_partyId < indexForDecreasingSize) {
        myNumOfElementsToExpect = numOfElementsForParties + 1;
    }
    for(int i=0;i<N;i++){

        if(i<=2*T){
            //recBufsBytes[i].resize(myNumOfElementsToExpect * fieldByteSize);
            recBufsElements[i].resize(myNumOfElementsToExpect);
        }
        else{
            //recBufsBytes[i].resize(0);
            recBufsElements[i].resize(0);
        }

    }



    roundFunctionSyncElements(sendBufsElements, recBufsElements,20);


    xyMinusR.resize(myNumOfElementsToExpect);
    xyMinusRBytes.resize(myNumOfElementsToExpect*fieldByteSize);

    //reconstruct the shares that I am responsible of recieved from the other parties
    vector<FieldType> xyMinurAllShares(2*T+1);
    counter = 0;

    for (int k = 0;k < myNumOfElementsToExpect; k++)//go over only the logit gates
    {
        for (int i = 0; i < N; i++) {

            if(i<=2*T) {

                xyMinurAllShares[counter] = recBufsElements[i][k];//field->bytesToElement(recBufsBytes[i].data() + (k * fieldByteSize));
                counter++;
            }

        }

        xyMinusR[k] = interpolateForParty(xyMinurAllShares, m_partyId);
        counter = 0;

    }

    //field->elementVectorToByteVector(xyMinusR, xyMinusRBytes);

    //prepare the send buffers
    for(int i=0; i<N; i++){
        sendBufsElements[i] = xyMinusR;
    }


    for(int i=0; i<N; i++){

        currentNumOfElements = numOfElementsForParties;
        if(i<indexForDecreasingSize)
            currentNumOfElements++;

        recBufsElements[i].resize(currentNumOfElements);

    }

    roundFunctionSyncElements(sendBufsElements, recBufsElements,21);


    xyMinusR.resize(acctualNumOfMultGates);
    counter = 0;

    for(int i=0; i<N; i++){

        currentNumOfElements = numOfElementsForParties;
        if(i<indexForDecreasingSize)
            currentNumOfElements++;

        //fill the send buf according to the number of elements to send to each party
        for(int j=0; j<currentNumOfElements; j++) {

            xyMinusR[counter] = recBufsElements[i][j];//field->bytesToElement(recBufsBytes[i].data() + (j * fieldByteSize));
            counter++;

        }

    }


    indexInRandomArray -= index;
    index = 0;

    //after the xPlusAAndYPlusB array is filled, we are ready to fill the output of the mult gates
    for (int k = circuit.getLayers()[currentCirciutLayer];
         k < circuit.getLayers()[currentCirciutLayer + 1]; k++)//go over only the logit gates
    {
        auto gate = circuit.getGates()[k];

        if (gate.gateType == MULT) {

            gateShareArr[gate.output] = randomTAnd2TShares[2*indexInRandomArray] + xyMinusR[index];

            index++;
            indexInRandomArray++;

        }
    }

    return index;
}
template <class FieldType>
int ProtocolParty<FieldType>::processMultDN(int indexInRandomArray) {

    int index = 0;
    int fieldByteSize = field->getElementSizeInBytes();
    int maxNumberOfLayerMult = circuit.getLayers()[currentCirciutLayer + 1] - circuit.getLayers()[currentCirciutLayer];
    vector<FieldType> xyMinusRShares(maxNumberOfLayerMult*2);//hold both in the same vector to send in one batch

    vector<FieldType> xyMinusR;//hold both in the same vector to send in one batch
    vector<byte> xyMinusRBytes;

    vector<vector<FieldType>> sendBufsElements(N);
    vector<vector<FieldType>> recBufsElements(N);



    //generate the shares for x+a and y+b. do it in the same array to send once
    for (int k = circuit.getLayers()[currentCirciutLayer];
         k < circuit.getLayers()[currentCirciutLayer + 1]; k++)//go over only the logit gates
    {
        auto gate = circuit.getGates()[k];

        if (gate.gateType == MULT) {

            //compute the share of xy-r
            xyMinusRShares[index] = gateShareArr[gate.input1]*gateShareArr[gate.input2] - randomTAnd2TShares[2*indexInRandomArray+1];


            indexInRandomArray++;

            index++;
        }
    }

    if(index==0)
        return 0;

    //set the acctual number of mult gate proccessed in this layer
    int acctualNumOfMultGates = index;
    int numOfElementsForParties = acctualNumOfMultGates/N;
    int indexForDecreasingSize = acctualNumOfMultGates - numOfElementsForParties *N;

    int counter=0;
    int currentNumOfElements;
    for(int i=0; i<N; i++){

        currentNumOfElements = numOfElementsForParties;
        if(i<indexForDecreasingSize)
            currentNumOfElements++;

        //fill the send buf according to the number of elements to send to each party
        sendBufsElements[i].resize(currentNumOfElements);
        for(int j=0; j<currentNumOfElements; j++) {

            sendBufsElements[i][j] = xyMinusRShares[counter];
            counter++;

        }

    }

    //resize the recbuf array.
    int myNumOfElementsToExpect = numOfElementsForParties;
    if (m_partyId < indexForDecreasingSize) {
        myNumOfElementsToExpect = numOfElementsForParties + 1;
    }
    for(int i=0;i<N;i++){

        //recBufsBytes[i].resize(myNumOfElementsToExpect*fieldByteSize);
        recBufsElements[i].resize(myNumOfElementsToExpect);
    }



    roundFunctionSyncElements(sendBufsElements, recBufsElements,20);


    xyMinusR.resize(myNumOfElementsToExpect);
    //xyMinusRBytes.resize(myNumOfElementsToExpect*fieldByteSize);

    //reconstruct the shares that I am responsible of recieved from the other parties
    vector<FieldType> xyMinurAllShares(N);

    for (int k = 0;k < myNumOfElementsToExpect; k++)//go over only the logit gates
    {
        for (int i = 0; i < N; i++) {

            xyMinurAllShares[i] = recBufsElements[i][k];////field->bytesToElement(recBufsBytes[i].data() + (k * fieldByteSize));
        }

        // reconstruct the shares by P0
        xyMinusR[k] = interpolate(xyMinurAllShares);
    }


    //prepare the send buffers
    for(int i=0; i<N; i++){

        sendBufsElements[i] = xyMinusR;
    }


    for(int i=0; i<N; i++){

        currentNumOfElements = numOfElementsForParties;
        if(i<indexForDecreasingSize)
            currentNumOfElements++;

        recBufsElements[i].resize(currentNumOfElements);

    }

    roundFunctionSyncElements(sendBufsElements, recBufsElements,21);


    xyMinusR.resize(acctualNumOfMultGates);
    counter = 0;

    for(int i=0; i<N; i++){

        currentNumOfElements = numOfElementsForParties;
        if(i<indexForDecreasingSize)
            currentNumOfElements++;

        //fill the send buf according to the number of elements to send to each party
        for(int j=0; j<currentNumOfElements; j++) {

            xyMinusR[counter] = recBufsElements[i][j];//field->bytesToElement(recBufsBytes[i].data() + (j * fieldByteSize));
            counter++;

        }

    }


    indexInRandomArray -= index;
    index = 0;

    //after the xPlusAAndYPlusB array is filled, we are ready to fill the output of the mult gates
    for (int k = circuit.getLayers()[currentCirciutLayer];
         k < circuit.getLayers()[currentCirciutLayer + 1]; k++)//go over only the logit gates
    {
        auto gate = circuit.getGates()[k];

        if (gate.gateType == MULT) {

            gateShareArr[gate.output] = randomTAnd2TShares[2*indexInRandomArray] + xyMinusR[index];

            index++;
            indexInRandomArray++;

        }
    }

    return index;
}


template <class FieldType>
void ProtocolParty<FieldType>::offlineDNForMultiplication(int numOfTriples){

    if(isPRF)
        alternativeGenerateRandom2TAndTShares(numOfTriples,randomTAnd2TShares);
    else
        generateRandom2TAndTShares(numOfTriples,randomTAnd2TShares);

}

template <class FieldType>
void ProtocolParty<FieldType>::verificationPhase() {

    //calc the number of times we need to run the verification -- ceiling
    int delta =   (5 + field->getElementSizeInBytes() - 1) / field->getElementSizeInBytes();


    int numOfOutputGates = circuit.getNrOfOutputGates();
    int numOfInputGates = circuit.getNrOfInputGates();
    int numOfMultGates = circuit.getNrOfMultiplicationGates();
    auto t1 = high_resolution_clock::now();

    //FieldType *multGatesShares = new FieldType[numOfMultGates*3];

    //1. prepare random values by generating a key
    //first generate the common aes key
    auto key = generateCommonKey();
    auto t2 = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(t2-t1).count();
//    if(flag_print_timings) {
//        cout << "time in milliseconds for generating key: " << duration << endl;
//    }


    t1 = high_resolution_clock::now();
    generatePseudoRandomElements(key, randomElementsForVerification, randomElementsForVerification.size());

    t2 = high_resolution_clock::now();

    duration = duration_cast<milliseconds>(t2-t1).count();
//    if(flag_print_timings) {
//        cout << "time in milliseconds for creating randoms: " << duration << endl;
//    }

    //1 (b) and 1(C)
    //called in the offline phase

    t1 = high_resolution_clock::now();
    int index = 0;
    for (int k = 0; k < numOfInputGates; k++) {

        auto gate = circuit.getGates()[k];

        if (gate.gateType == INPUT) {
            wSharesForVerification[index] = gateShareArr[gate.output];

            index++;
        }
    }

    index = 0;
    for (int k = numOfInputGates - 1; k < M - numOfOutputGates + 1; k++) {

        auto gate = circuit.getGates()[k];
        if (gate.gateType == MULT) {

            multGatesSharesForVerification[2*index] = gateShareArr[gate.input1]*gateShareArr[gate.input2];

            //multGatesShares[2*index+1] = gateShareArr[gate.input2];

            multGatesSharesForVerification[2*index+1] = gateShareArr[gate.output];
            index++;
        }


    }

    t2 = high_resolution_clock::now();

    duration = duration_cast<milliseconds>(t2-t1).count();
//    if(flag_print_timings) {
//        cout << "time in milliseconds for preparing arrays: " << duration << endl;
//    }


    //verify that all sharings on all wires are of degree t
    //2 (a)
    vector<FieldType> uShares(delta);

    index = 0;

    for(int iter=0; iter<delta; iter++) {
        uShares[iter] = *field->GetZero();
        for (int i = 0; i < numOfInputGates; i++) {

            uShares[iter] += randomElementsForVerification[index] * wSharesForVerification[i];
            index++;

        }


        for (int i = 0; i < numOfMultGates; i++) {

            uShares[iter] += randomElementsForVerification[index]*multGatesSharesForVerification[2*i + 1];
        }

        uShares[iter] += secureRandomTShares[randomOffset];
        randomOffset++;
    }

    vector<FieldType> uSecrets(delta);
    t1 = high_resolution_clock::now();
    openShare(delta,uShares, uSecrets, T);
    t2 = high_resolution_clock::now();

    duration = duration_cast<milliseconds>(t2-t1).count();
//    if(flag_print_timings) {
//        cout << "time in milliseconds for open share: " << duration << endl;
//    }

    //verify that sll multiplications were valid
    //3 (a)
    vector<FieldType> v2TShares(delta);
    vector<FieldType> vSecrets(delta);


    t1 = high_resolution_clock::now();
    for(int iter=0; iter<delta; iter++) {
        v2TShares[iter] = *field->GetZero();
        for (int i = 0; i < numOfMultGates; i++) {

            v2TShares[iter] += randomElementsForVerification[index] * (multGatesSharesForVerification[i*2] - multGatesSharesForVerification[i*2+1]);
            index++;

        }



        v2TShares[iter] += secureRandom2TShares[randomOffset];
        randomOffset++;
    }
    t2 = high_resolution_clock::now();

    duration = duration_cast<milliseconds>(t2-t1).count();
//    if(flag_print_timings) {
//        cout << "time in milliseconds for mults: " << duration << endl;
//    }

    //3 (b)
    openShare(delta, v2TShares, vSecrets, 2*T);

    //3 (c)
    vector<FieldType> w(delta);
    randomOffset = randomOffset - delta;
    for(int i=0; i<delta; i++){

        w[i] = secureRandomTShares[randomOffset] - vSecrets[i];
        randomOffset++;
    }

    //3 (d)
    if (checkZero(delta, w)==false)
        exit(0);

  }

template <class FieldType>
bool ProtocolParty<FieldType>::checkZero(int numOfChecks, vector<FieldType> &toCheck) {

    vector<FieldType> rv2TShares(numOfChecks);
    vector<FieldType> rvSecrets(numOfChecks);
    // 2. compute locally [r*v]_2T = [r]_T * [v]_T

    for(int i=0; i<numOfChecks; i++){

        //take the t share only
        rv2TShares[i] = secureRandomTShares[randomOffset]*toCheck[i];

        //increment the pointer and "waist" the corresponding 2T random share
        randomOffset++;
    }

    //3. open the share
    openShare(numOfChecks, rv2TShares, rvSecrets, 2*T);

    //4 check equality to zero

    for(int i=0; i<numOfChecks; i++){

        if(rvSecrets[i] != *field->GetZero()) {
            cout<<"bassssssaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa not equal to zero"<<endl;
            return false;
        }
    }

    return true;


}

  template <class FieldType>
  vector<byte> ProtocolParty<FieldType>::generateCommonKey(){

      int fieldByteSize = field->getElementSizeInBytes();

      //calc the number of elements needed for 128 bit AES key
      int numOfRandomShares = 16/field->getElementSizeInBytes() + 1;
      //vector<FieldType> randomSharesArray(numOfRandomShares);
      vector<FieldType> aesArray(numOfRandomShares);
      vector<byte> aesKey(numOfRandomShares*fieldByteSize);


      //get the secure random shares
      vector<FieldType> secRandShares(secureRandomTShares.begin() + randomOffset,secureRandomTShares.begin() + randomOffset + numOfRandomShares);
      randomOffset+=numOfRandomShares;



      openShare(numOfRandomShares, secRandShares, aesArray, T);


      //turn the aes array into bytes to get the common aes key.
      for(int i=0; i<numOfRandomShares;i++){

          for(int j=0; j<numOfRandomShares;j++) {
              field->elementToBytes(aesKey.data() + (j * fieldByteSize), aesArray[j]);
          }
      }

      //reduce the size of the key to 16 bytes
      aesKey.resize(16);

      return aesKey;

  }

  template <class FieldType>
  void ProtocolParty<FieldType>::openShare(int numOfRandomShares, vector<FieldType> &Shares, vector<FieldType> &secrets, int degree){

      vector<vector<FieldType>> sendBufsElements(N);
      vector<vector<FieldType>> recBufsElements(N);

      vector<FieldType> x1(N);
      int fieldByteSize = field->getElementSizeInBytes();

      //calc the number of elements needed for 128 bit AES key

      //resize vectors
      for(int i=0; i < N; i++)
      {
          sendBufsElements[i].resize(numOfRandomShares);
          recBufsElements[i].resize(numOfRandomShares);
      }


      //copy the same data for all parties
      for(int i=0; i<N; i++){

          sendBufsElements[i] = Shares;
      }

      //call the round function to send the shares to all the users and get the other parties share
      roundFunctionSyncElements(sendBufsElements, recBufsElements,12);

      //reconstruct each set of shares to get the secret

      for(int k=0; k<numOfRandomShares; k++){

          //get the set of shares for each element
          for(int i=0; i < N; i++) {

              x1[i] = recBufsElements[i][k];//field->bytesToElement(recBufsBytes[i].data() + (k*fieldByteSize));
          }

          secrets[k] = reconstructShare(x1, degree);
      }

  }


template <class FieldType>
void ProtocolParty<FieldType>::generatePseudoRandomElements(vector<byte> & aesKey, vector<FieldType> &randomElementsToFill, int numOfRandomElements){

    int fieldSize = field->getElementSizeInBytes();
    int fieldSizeBits = field->getElementSizeInBits();
    bool isLongRandoms;
    int size;

    if(fieldSize==1){
        size = 1;
    }

    else if(fieldSize>4){
      isLongRandoms = true;
      size = 8;
    }
    else{

      isLongRandoms = false;
      size = 4;
    }

    if (flag_print) {
        cout << "size is" << size << "for party : " << m_partyId;
    }

    auto t1 = high_resolution_clock::now();
    PrgFromOpenSSLAES prg((numOfRandomElements*size/16) + 1);
    SecretKey sk(aesKey, "aes");
    prg.setKey(sk);


    byte *randBytes;
    prg.getPRGBytes(randBytes, numOfRandomElements * size);
    auto t2 = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(t2-t1).count();
    if(flag_print_timings) {
        cout << "time in milliseconds for getting random bytes: " << duration << endl;
    }

    t1 = high_resolution_clock::now();

    if(fieldSize==1){
        //GF2_8, just copy the bytes to the elements array
        memcpy(randomElementsToFill.data() , randBytes, numOfRandomElements * size);
    }

    else if(isLongRandoms) {
        //mersenne61 shift left before copying in order the make the bytes fit the field (there is
        // one exception where the bytes might ad up to be p, but since we mult and adjust, the result will
        // be back in the field
        for(int i=0; i<numOfRandomElements; i++){
            randBytes[8*i+7] = randBytes[8*i+7]>>(64 - fieldSizeBits);
        }
        //copy to the element array after ajdustment
        memcpy(randomElementsToFill.data() , randBytes, numOfRandomElements * size);

    }
      else{
        //mersenne31 shift left before copying in order the make the bytes fit the field (there is
        // one exception where the bytes might ad up to be p, but since we mult and adjust, the result will
        // be back in the field
        for(int i=0; i<numOfRandomElements; i++){
            randBytes[4*i+3] = randBytes[4*i+3]>>(32 - fieldSizeBits);
        }
        //copy to the element array after ajdustment
        memcpy(randomElementsToFill.data() , randBytes, numOfRandomElements * size);
    }

    t2 = high_resolution_clock::now();
    duration = duration_cast<milliseconds>(t2-t1).count();
//    if(flag_print_timings) {
//        cout << "time in milliseconds for setting to elements: " << duration << endl;
//    }
}


/**
 * the function Walk through the circuit and reconstruct output gates.
 * @param circuit
 * @param gateShareArr
 * @param alpha
 */
template <class FieldType>
void ProtocolParty<FieldType>::outputPhase()
{
    int count=0;
    vector<FieldType> x1(N); // vector for the shares of my outputs
    vector<vector<FieldType>> sendBufsElements(N);
    vector<vector<FieldType>> recBufElements(N);

    FieldType num;
    ofstream myfile;
    myfile.open(outputFile);

    for(int k=M-numOfOutputGates; k < M; k++)
    {
        if(circuit.getGates()[k].gateType == OUTPUT)
        {
            // send to party (which need this gate) your share for this gate
            sendBufsElements[circuit.getGates()[k].party].push_back(gateShareArr[circuit.getGates()[k].input1]);
        }
    }


    int fieldByteSize = field->getElementSizeInBytes();
    for(int i=0; i < N; i++)
    {
        recBufElements[i].resize(sendBufsElements[m_partyId].size());
    }

    roundFunctionSyncElements(sendBufsElements, recBufElements,7);



    int counter = 0;
    if(flag_print) {
        cout << "endnend" << endl;}
    for(int k=M-numOfOutputGates ; k < M; k++) {
        if(circuit.getGates()[k].gateType == OUTPUT && circuit.getGates()[k].party == m_partyId)
        {
            for(int i=0; i < N; i++) {

                x1[i] = recBufElements[i][counter];//field->bytesToElement(recBufBytes[i].data() + (counter*fieldByteSize));
            }


            // my output: reconstruct received shares
            if (!checkConsistency(x1, T))
            {
                // someone cheated!
                //if(flag_print) {
                    cout << "cheating!!!" << '\n';//}
                return;
            }
            if(flag_print_output)
                cout << "the result for "<< circuit.getGates()[k].input1 << " is : " << field->elementToString(interpolate(x1)) << '\n';


            counter++;
        }
    }

    // close output file
    myfile.close();
}


template <class FieldType>
void ProtocolParty<FieldType>::roundFunctionSyncElements(vector<vector<FieldType>> &sendBufs, vector<vector<FieldType>> &recBufs, int round) {

    //cout<<"in roundFunctionSync "<< round<< endl;

    int numThreads = 10;//parties.size();
    int numPartiesForEachThread;

    if (parties.size() <= numThreads){
        numThreads = parties.size();
        numPartiesForEachThread = 1;
    } else{
        numPartiesForEachThread = (parties.size() + numThreads - 1)/ numThreads;
    }


    recBufs[m_partyId] = move(sendBufs[m_partyId]);
    //recieve the data using threads
    vector<thread> threads(numThreads);
    for (int t=0; t<numThreads; t++) {
        if ((t + 1) * numPartiesForEachThread <= parties.size()) {
            threads[t] = thread(&ProtocolParty::exchangeDataElements, this, ref(sendBufs), ref(recBufs),
                                t * numPartiesForEachThread, (t + 1) * numPartiesForEachThread);
        } else {
            threads[t] = thread(&ProtocolParty::exchangeDataElements, this, ref(sendBufs), ref(recBufs), t * numPartiesForEachThread, parties.size());
        }
    }
    for (int t=0; t<numThreads; t++){
        threads[t].join();
    }

}


template <class FieldType>
void ProtocolParty<FieldType>::exchangeDataElements(vector<vector<FieldType>> &sendBufs, vector<vector<FieldType>> &recBufs, int first, int last) {


    //cout<<"in exchangeData";
    for (int i = first; i < last; i++) {

        if ((m_partyId) < parties[i]->getID()) {


            if (sendBufs[parties[i]->getID()].size() > 0) {
                //send shares to my input bits
                parties[i]->getChannel()->write((byte *) sendBufs[parties[i]->getID()].data(),
                                                sendBufs[parties[i]->getID()].size() * field->getElementSizeInBytes());
                //cout<<"write the data:: my Id = " << m_partyId - 1<< "other ID = "<< parties[i]->getID() <<endl;
            }

            if (recBufs[parties[i]->getID()].size() > 0) {
                //receive shares from the other party and set them in the shares array
                parties[i]->getChannel()->read((byte *) recBufs[parties[i]->getID()].data(),
                                               recBufs[parties[i]->getID()].size() * field->getElementSizeInBytes());
                //cout<<"read the data:: my Id = " << m_partyId-1<< "other ID = "<< parties[i]->getID()<<endl;
            }

        } else {

            if (recBufs[parties[i]->getID()].size() > 0) {
                //receive shares from the other party and set them in the shares array
                parties[i]->getChannel()->read((byte *) recBufs[parties[i]->getID()].data(),
                                               recBufs[parties[i]->getID()].size() * field->getElementSizeInBytes());
                //cout<<"read the data:: my Id = " << m_partyId-1<< "other ID = "<< parties[i]->getID()<<endl;
            }

            if (sendBufs[parties[i]->getID()].size() > 0) {

                //send shares to my input bits
                parties[i]->getChannel()->write((byte *) sendBufs[parties[i]->getID()].data(),
                                                sendBufs[parties[i]->getID()].size() * field->getElementSizeInBytes());
                //cout<<"write the data:: my Id = " << m_partyId-1<< "other ID = "<< parties[i]->getID() <<endl;
            }

        }

    }
}







template <class FieldType>
void ProtocolParty<FieldType>::roundFunctionSync(vector<vector<byte>> &sendBufs, vector<vector<byte>> &recBufs, int round) {

    //cout<<"in roundFunctionSync "<< round<< endl;

    int numThreads = 10;//parties.size();
    int numPartiesForEachThread;

    if (parties.size() <= numThreads){
        numThreads = parties.size();
        numPartiesForEachThread = 1;
    } else{
        numPartiesForEachThread = (parties.size() + numThreads - 1)/ numThreads;
    }


    recBufs[m_partyId] = move(sendBufs[m_partyId]);
    //recieve the data using threads
    vector<thread> threads(numThreads);
    for (int t=0; t<numThreads; t++) {
        if ((t + 1) * numPartiesForEachThread <= parties.size()) {
            threads[t] = thread(&ProtocolParty::exchangeData, this, ref(sendBufs), ref(recBufs),
                                t * numPartiesForEachThread, (t + 1) * numPartiesForEachThread);
        } else {
            threads[t] = thread(&ProtocolParty::exchangeData, this, ref(sendBufs), ref(recBufs), t * numPartiesForEachThread, parties.size());
        }
    }
    for (int t=0; t<numThreads; t++){
        threads[t].join();
    }

}


template <class FieldType>
void ProtocolParty<FieldType>::exchangeData(vector<vector<byte>> &sendBufs, vector<vector<byte>> &recBufs, int first, int last){


    //cout<<"in exchangeData";
    for (int i=first; i < last; i++) {

        if ((m_partyId) < parties[i]->getID()) {


            //send shares to my input bits
            parties[i]->getChannel()->write(sendBufs[parties[i]->getID()].data(), sendBufs[parties[i]->getID()].size());
            //cout<<"write the data:: my Id = " << m_partyId - 1<< "other ID = "<< parties[i]->getID() <<endl;
            
            //receive shares from the other party and set them in the shares array
            parties[i]->getChannel()->read(recBufs[parties[i]->getID()].data(), recBufs[parties[i]->getID()].size());
            //cout<<"read the data:: my Id = " << m_partyId-1<< "other ID = "<< parties[i]->getID()<<endl;

        } else{


            //receive shares from the other party and set them in the shares array
            parties[i]->getChannel()->read(recBufs[parties[i]->getID()].data(), recBufs[parties[i]->getID()].size());
            //cout<<"read the data:: my Id = " << m_partyId-1<< "other ID = "<< parties[i]->getID()<<endl;

            //send shares to my input bits
            parties[i]->getChannel()->write(sendBufs[parties[i]->getID()].data(), sendBufs[parties[i]->getID()].size());
            //cout<<"write the data:: my Id = " << m_partyId-1<< "other ID = "<< parties[i]->getID() <<endl;

        }
    }
}


template <class FieldType>
ProtocolParty<FieldType>::~ProtocolParty()
{
    protocolTimer->writeToFile();
    delete protocolTimer;
    delete field;
    delete timer;
}


#endif /* PROTOCOLPARTY_H_ */
