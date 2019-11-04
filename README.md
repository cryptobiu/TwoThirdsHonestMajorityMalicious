### Two-Thirds Honest-Majority MPC for Malicious Adversaries at Almost the Cost of Semi-Honest

The repository implements the Two-Thirds Honest-Majority MPC for Malicious Adversaries at Almost the Cost of Semi-Honest [article](https://eprint.iacr.org/2019/658.pdf).  

##### Abstract
Abstract. Secure multiparty computation (MPC) enables a set of parties to securely carry out a joint computation of their private inputs
without revealing anything but the output. Protocols for semi-honest
adversaries guarantee security as long as the corrupted parties run the
specified protocol and ensure that nothing is leaked in the transcript.
In contrast, protocols for malicious adversaries guarantee security in the
presence of arbitrary adversaries who can run any attack strategy. Security for malicious adversaries is typically what is needed in practice (and
is always preferred), but comes at a significant cost.
In this paper, we present the first protocol for a two-thirds honest majority that achieves security in the presence of malicious adversaries at essentially the exact same cost as the best known protocols for semi-honest
adversaries. Our construction is not a general transformation and thus
it is possible that better semi-honest protocols will be constructed which
do not support our transformation. Nevertheless, for the current state-ofthe-art for many parties (based on Shamir sharing), our protocol invokes
the best semi-honest multiplication protocol exactly once per multiplication gate (plus some additional local computation that is negligible to
the overall cost). Concretely, the best version of our protocol requires
each party to send on average of just 2 2
3 elements per multiplication
gate (when the number of multiplication gates is at least the number
of parties). This is four times faster than the previous-best protocol of
Barak et al. (ACM CCS 2018) for small fields, and twice as fast as the
previous-best protocol of Chida et al. (CRYPTO 2018) for large fields.

##### Note

The code here is close to the specified protocol, but there are slight differences that have security ramifications.
Thus, this code can be used for efficiency comparisons to other protocols, but is not secure as is. 


##### Installation

The protocol written in c++ and uses c++11 standard. It uses [libscapi](https://github.com/cryptobiu/libscapi).  
For `libscapi` installation instructions, visit [here](https://github.com/cryptobiu/libscapi/blob/master/build_scripts/INSTALL.md).  
After you installed `libscapi`, run `cmake . && make`

##### Usage

The protocol designed for at least 3 parties.
To run the the protocol open a terminal and run:  
`run_protocol.sh <min_party_id> <max_party_id> <number_of_parties> <input_file> <circuit_file> <filed type> <parties_file> <-is_honest> <number_of_iterations> <is_prf>` 

* field_type can be one of this values:
    * ZpMersenne31
    * ZpMersenne61
    * GF2_8LookupTable

* parties_file - a file that contains the ip addresses and the port of all the parties. An example file can be found [here](../master/Parties.txt).


    