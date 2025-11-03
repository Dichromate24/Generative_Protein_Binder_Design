import requests
import json
import time
import os
from pathlib import Path

import py3Dmol
import webbrowser
import tempfile
import zipfile
import io


def AF2_NIM(API_KEY, sequence):
    AF2_endpoint = "https://health.api.nvidia.com/v1/biology/deepmind/alphafold2"
    AF2_status_endpoint = "https://health.api.nvidia.com/v1/status" 

    AF2_header = {
        "Authorization": f"Bearer {API_KEY}",
        "Content-Type": "application/json",
        "Accept": "application/json",
        #"NVCF-POLL-SECONDS": "600"
    }

    AF2_payload = {
        "sequence": sequence,
        "algorithm": "mmseqs2"
    }

    AF2_response = requests.post(
                    AF2_endpoint,
                    headers=AF2_header,
                    json=AF2_payload,
                    timeout=(10,310)
                )

    if AF2_response.status_code == 200:    
        AF2_result = AF2_response.json()
        return AF2_result
    
    elif AF2_response.status_code == 202:
        print("Polling...")

        req_id = AF2_response.headers.get("NVCF-REQID") or AF2_response.headers.get("nvcf-reqid")
        if not req_id:
            return f"Error: 202 without NVCF-REQID header. Headers: {dict(AF2_response.headers)}"
        print(f"request id: {req_id}")

        max_attempt = 180
        poll_sec = 10
        for attempt in range(max_attempt):
            print(f"trying attempt {attempt + 1}...")
            
            try:
                time.sleep(poll_sec)
                result = requests.get(f"{AF2_status_endpoint}/{req_id}", 
                                      headers=AF2_header, 
                                      timeout=(10,310))

                if result.status_code == 200:
                    AF2_result = result.json()
                    return AF2_result
                elif result.status_code == 202:
                    continue
                else:
                    return f"Status polling failed: {result.status_code} - {result.text}"
                
            except requests.exceptions.ReadTimeout:
                print("Read timed out (long-poll). Continuing to poll...")
                time.sleep(poll_sec)
                continue
            except requests.exceptions.RequestException as e:
                # Network hiccup; keep trying unless near the end
                print(f"Request exception while polling: {e}")
                time.sleep(poll_sec)
                continue
        return f"Timeout after {max_attempt * poll_sec // 60} minutes waiting for result (req_id={req_id})."

    else:
            print(f"error status code: {AF2_response.status_code}")
            return f"Error: {AF2_response.status_code} - {AF2_response.text}"


def OF3_NIM(API_KEY, example_input, protein_name, msa):

    OF3_endpoint = "https://health.api.nvidia.com/v1/biology/openfold/openfold3/predict"

    OF3_header = {
    "content-type": "application/json",
    "Authorization": f"Bearer {API_KEY}",
    "NVCF-POLL-SECONDS": "300",
    }

    OF3_payload = {
        "request_id": protein_name,
        "inputs": [
            {
                "input_id": protein_name,
                "molecules": [
                    {
                        "type": "protein",
                        "id": "A",
                        "sequence": f"{example_input}",
                        "msa": {
                            "main_db": {
                                "csv": {
                                    "alignment": msa,
                                    "format": "csv",
                                }
                            },
                        },
                    },
                ],
                "output_format": "pdb",
            }
        ],
    }

    OF3_response = requests.post(OF3_endpoint, headers=OF3_header, json=OF3_payload)

    if OF3_response.status_code == 200:
        return OF3_response.text
    else:
        print(f"Unexpected HTTP status: {OF3_response.status_code}")
        print(f"Response: {OF3_response.text}")


def ESM_NIM(API_KEY, sequence):
    ESM_endpoint = "https://health.api.nvidia.com/v1/biology/nvidia/esmfold"

    ESM_header = {
        "Authorization": f"Bearer {API_KEY}",
        "Content-Type": "application/json",
        "Accept": "application/json",
    }

    ESM_payload = {
        "sequence": sequence
    }

    ESM_response = requests.post(
        ESM_endpoint,
        headers=ESM_header,
        json=ESM_payload
    )

    return ESM_response.json()



def RFD_NIM(API_KEY, sequence):
    RFD_endpoint = "https://health.api.nvidia.com/v1/biology/ipd/rfdiffusion/generate"

    RFD_header = {
        "Authorization": f"Bearer {API_KEY}",
        "Content-Type": "application/json",
        "Accept": "application/json",
        #"NVCF-POLL-SECONDS": "600"
    }

    RFD_payload = {
        "input_pdb": sequence,
        "contigs": "A1-25/0 70-100",
        "hotspot_res": ["A14","A15","A17","A18"],
        "diffusion_steps": 15,
    }

    RFD_response = requests.post(
                        RFD_endpoint,
                        headers = RFD_header,
                        json = RFD_payload
                    )
    
    return RFD_response.text


def PMPNN_NIM(API_KEY, sequence):
    PMPNN_endpoint = "https://health.api.nvidia.com/v1/biology/ipd/proteinmpnn/predict"

    PMPNN_header = {
        "Authorization": f"Bearer {API_KEY}",
        "Content-Type": "application/json",
        "Accept": "application/json",
        #"NVCF-POLL-SECONDS": "600"
    }

    PMPNN_payload = {
        "input_pdb": sequence,
        "input_pdb_chains" : ["A"],
        "ca_only" : False,
        "use_soluble_model" : False,
        "num_seq_per_target" : 10,
        "sampling_temp" : [0.1]
    }

    PMPNN_response = requests.post(
                        PMPNN_endpoint,
                        headers = PMPNN_header,
                        json = PMPNN_payload
                    )
    
    return PMPNN_response.text


def MULTIMER_NIM(API_KEY, binder_pair):

    MULTIMER_endpoint = "https://health.api.nvidia.com/v1/biology/deepmind/alphafold2-multimer"
    MULTIMER_status_endpoint = "https://health.api.nvidia.com/v1/status"

    MULTIMER_header = {
        "Authorization": f"Bearer {API_KEY}",
        "Content-Type": "application/json",
        "Accept": "application/json",
        #"NVCF-POLL-SECONDS": "600"
    }

    MULTIMER_payload = {
        "sequences" : binder_pair,
        "selected_models" : [1],
        "relax_prediction": False,         # skip relaxation, faster
        "databases": ["small_bfd"]         # small dataase

        #"algorithm": "jackhmmer",
        #"e_value": 0.0001,
        #"iterations": 1,
        #"databases": ["uniref90", "small_bfd", "mgnify"],
        #"relax_prediction": True,
    }

    MULTIMER_response = requests.post(
                    MULTIMER_endpoint,
                    headers=MULTIMER_header,
                    json=MULTIMER_payload,
                    timeout=(10,310)
                )
    
    if MULTIMER_response.status_code == 200:    
        return MULTIMER_response
    
    elif MULTIMER_response.status_code == 202:
        print("Polling...")

        req_id = MULTIMER_response.headers.get("NVCF-REQID") or MULTIMER_response.headers.get("nvcf-reqid")
        if not req_id:
            return f"Error: 202 without NVCF-REQID header. Headers: {dict(MULTIMER_response.headers)}"
        print(f"request id: {req_id}")

        max_attempt = 180
        poll_sec = 10
        for attempt in range(max_attempt):
            print(f"trying attempt {attempt + 1}...")
            
            try:
                time.sleep(poll_sec)
                result = requests.get(f"{MULTIMER_status_endpoint}/{req_id}", 
                                      headers=MULTIMER_header, 
                                      timeout=(10,310))

                if result.status_code == 200:
                    return result
                elif result.status_code == 202:
                    continue
                else:
                    return f"Status polling failed: {result.status_code} - {result.text}"
                
            except requests.exceptions.ReadTimeout:
                print("Read timed out (long-poll). Continuing to poll...")
                time.sleep(poll_sec)
                continue
            except requests.exceptions.RequestException as e:
                # Network hiccup; keep trying unless near the end
                print(f"Request exception while polling: {e}")
                time.sleep(poll_sec)
                continue
        return f"Timeout after {max_attempt * poll_sec // 60} minutes waiting for result (req_id={req_id})."

    else:
            print(f"error status code: {MULTIMER_response.status_code}")
            return f"Error: {MULTIMER_response.status_code} - {MULTIMER_response.text}"







def visualize_pdb_py3dmol(pdb_text: str, width: int = 800, height: int = 600, spin: bool = False) -> str:
    """
    Make an HTML file with a 3Dmol viewer and open it in the browser.
    Returns the HTML file path.
    """
    view = py3Dmol.view(width=width, height=height)
    view.addModel(pdb_text, 'pdb')

    view.setStyle({'cartoon': {'color': 'spectrum'}})  # try 'chain' or 'ss' too

    view.zoomTo()
    if spin:
        view.spin(True)

    html = view._make_html()
    html_path = os.path.join(tempfile.gettempdir(), "pdb_viewer.html")
    with open(html_path, "w", encoding="utf-8") as f:
        f.write(html)

    webbrowser.open(f"file://{html_path}")
    return html_path



# Function to calculate average pLDDT over all residues 
def calculate_average_pLDDT(pdb_string):
    total_pLDDT = 0.0
    atom_count = 0
    pdb_lines = pdb_string.splitlines()
    for line in pdb_lines:
        # PDB atom records start with "ATOM"
        if line.startswith("ATOM"):
            atom_name = line[12:16].strip() # Extract atom name
            if atom_name == "CA":  # Only consider atoms with name "CA"
                try:
                    # Extract the B-factor value from columns 61-66 (following PDB format specifications)
                    pLDDT = float(line[60:66].strip())
                    total_pLDDT += pLDDT
                    atom_count += 1
                except ValueError:
                    pass  # Skip lines where B-factor can't be parsed as a float

    if atom_count == 0:
        return 0.0  # Return 0 if no N atoms were found

    average_pLDDT = total_pLDDT / atom_count
    return average_pLDDT





def main():

    print("\n\n")
    print("***************************************************************")
    print("************* GENERATIVE PROTEIN BINDER DESIGN ****************")
    print("***************************************************************")
    print("\n\n")

    #API_KEY = input("Please enter your API_KEY: ")
    #API_KEY = "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"

    # "Insulin B-chain (30 AA)": "FVNQHLCGSHLVEALYLVCGERGFFYTPKT"
    #example_input = "FVNQHLCGSHLVEALYLVCGERGFFYTPKT"
    #protein_name = "Insulin_B_chain"

    # ACE2 Receptor
    #example_input="STIEEQAKTFLDKFNHEAEDLFYQSSLASWNYNTNITEENVQNMNNAGDKWSAFLKEQSTLAQMYPLQEIQNLTVKLQLQALQQNGSSVLSEDKSKRLNTILNTMSTIYSTGKVCNPDNPQECLLLEPGLNEIMANSLDYNERLWAWESWRSEVGKQLRPLYEEYVVLKNEMARANHYEDYGDYWRGDYEVNGVDGYDYSRGQLIEDVEHTFEEIKPLYEHLHAYVRAKLMNAYPSYISPIGCLPAHLLGDMWGRFWTNLYSLTVPFGQKPNIDVTDAMVDQAWDAQRIFKEAEKFFVSVGLPNMTQGFWENSMLTDPGNVQKAVCHPTAWDLGKGDFRILMCTKVTMDDFLTAHHEMGHIQYDMAYAAQPFLLRNGANEGFHEAVGEIMSLSAATPKHLKSIGLLSPDFQEDNETEINFLLKQALTIVGTLPFTYMLEKWRWMVFKGEIPKDQWMKKWWEMKREIVGVVEPVPHDETYCDPASLFHVSNDYSFIRYYTRTLYQFQFQEALCQAAKHEGPLHKCDISNSTEAGQKLFNMLRLGKSEPWTLALENVVGAKNMNVRPLLNYFEPLFTWLKDQNKNSFVGWSTDWSPYAD"
    #protein_name = "ACE2_Receptor"

    # 5tpn
    #example_input = "NITEEFYQSTCSAVSKGYLSALRTGWYTSVITIELSNIKKIKCNGTDAKIKLIKQELDKYKNAVTELQLLMQSTPATNNQARGSGSGRSLGFLLGVGSAIASGVAVSKVLHLEGEVNKIKSALLSTNKAVVSLSNGVSVLTSKVLDLKNYIDKQLLPIVNKQSCSIPNIETVIEFQQKNNRLLEITREFSVNAGVTTPVSTYMLTNSELLSLINDMPITNDQKKLMSNNVQIVRQQSYSIMSIIKEEVLAYVVQLPLYGVIDTPCWKLHTSPLCTTNTKEGSNICLTRTDRGWYCDNAGSVSFFPQAETCKVQSNRVFCDTMNSLTLPSEVNLCNVDIFNPKYDCKIMTSKTDVSSSVITSLGAIVSCYGKTKCTASNKNRGIIKTFSNGCDYVSNKGVDTVSVGNTLYYVNKQEGKSLYVKGEPIINFYDPLVFPSDQFDASISQVNEKINQSLAFIRKSDELLSAIGGYIPEAPRDGQAYVRKDGEWVLLSTFLGGLVPRGSHHHHHH"
    #protein_name = "5tpn"
    '''

    NITEEFYQSTCSAVSKGYLSALRTGWYTSVITIELSNIKKIKCNGTDAKIKLIKQELDKYKNAVTELQLLMQSTPATNNQARGSGSGRSLGFLLGVGSAIASGVAVSKVLHLEGEVNKIKSALLSTNKAVVSLSNGVSVLTSKVLDLKNYIDKQLLPIVNKQSCSIPNIETVIEFQQKNNRLLEITREFSVNAGVTTPVSTYMLTNSELLSLINDMPITNDQKKLMSNNVQIVRQQSYSIMSIIKEEVLAYVVQLPLYGVIDTPCWKLHTSPLCTTNTKEGSNICLTRTDRGWYCDNAGSVSFFPQAETCKVQSNRVFCDTMNSLTLPSEVNLCNVDIFNPKYDCKIMTSKTDVSSSVITSLGAIVSCYGKTKCTASNKNRGIIKTFSNGCDYVSNKGVDTVSVGNTLYYVNKQEGKSLYVKGEPIINFYDPLVFPSDQFDASISQVNEKINQSLAFIRKSDELLSAIGGYIPEAPRDGQAYVRKDGEWVLLSTFLGGLVPRGSHHHHHH
    
    '''

    API_KEY = input("Please enter your API_KEY: ").strip()
    example_input = input("Please enter your protein sequence: ").strip()
    protein_name = input("Please enter your protein name: ").strip()
    model_name = input("Please enter your initial model (AF2/ OF3/ ESM): ").strip().upper()

    output_dir = Path(f"{protein_name}_{model_name}_output")
    output_dir.mkdir(parents=True, exist_ok=True)

    # ----------------  RUN AF2 NIM ---------------------------

    if model_name == "AF2":
        print("Running AlphaFold2...")
        AF2_result = AF2_NIM(API_KEY, example_input)
        
        # unpack AF2 results
        # result in list of 5 predicted structure, so concatenate with newlines
        all_pdb_text = "\n".join(x for x in AF2_result if isinstance(x, str))
        first_pdb = AF2_result[0]   
        
        # save all 5 AF2 results into pdb file
        (output_dir / f"{protein_name}_all_AF2_structures.pdb").write_text(all_pdb_text, encoding="utf-8")
        # save first AF2 result
        (output_dir / f"{protein_name}_first_structure.pdb").write_text(first_pdb, encoding="utf-8")

        # visualise AF2 results
        html_path = visualize_pdb_py3dmol(first_pdb)
        print("Opened viewer for AlphaFold2:", html_path)

    # ----------------  RUN OPENFOLD3 NIM -------------------
    elif model_name == "OF3":
        print("Running OpenFold3...")

        # openfold 3 requres MSA alignment in CSV format
        # add more lines if have more sequence for MSA
        # -1 is query sequence, 0,1,2... are for additional sequences
        msa_alignment_csv = (
            "key,sequence\n"
            f"-1,{example_input}"
        )

        OF3_result = OF3_NIM(API_KEY, example_input, protein_name, msa_alignment_csv)
        OF3_result = json.loads(OF3_result)

        pdb_text = (
            OF3_result.get("outputs", [{}])[0]
                .get("structures_with_scores", [{}])[0]
                .get("structure", "")
            )
        
        if not pdb_text.strip().startswith(("ATOM","HETATM","MODEL","HEADER")):
            raise ValueError("PDB text not found where expected.")

        filename = f"{protein_name}_first_structure.pdb"
        (output_dir / filename).write_text(pdb_text, encoding="utf-8")

        # visualise OpenFold 3 results
        html_path = visualize_pdb_py3dmol(pdb_text)
        print("Opened viewer for OpenFold3:", html_path)

    # ----------------  RUN ESMFOLD NIM -------------------
    elif model_name == "ESM":
        print("Running ESMFold...")

        ESM_result = ESM_NIM(API_KEY, example_input)
        pdb_text = "".join(ESM_result["pdbs"])

        filename = f"{protein_name}_first_structure.pdb"
        (output_dir / filename).write_text(pdb_text, encoding="utf-8")

        # visualise OpenFold 3 results
        html_path = visualize_pdb_py3dmol(pdb_text)
        print("Opened viewer for ESMFold:", html_path)

    else:
        raise ValueError("Model has to be 'AF2', 'OF3' or 'ESM'.")
    

    # ----------------  RUN RFDIFFUSION NIM -------------------
    print("Running RF-DIFFUSION...")

    # use the first output of AF2 OR OPENFOLD(either access via file, or directly)
    #RFD_input = first_AF2
    RFD_input = Path(f"{output_dir}/{protein_name}_first_structure.pdb")

    RFD_input = filter(lambda line: line.startswith("ATOM"), RFD_input.read_text().split("\n"))
    RFD_input = "\n".join(list(RFD_input)[:400])

    RFD_response = RFD_NIM(API_KEY, RFD_input)
    RFD_response = json.loads(RFD_response)["output_pdb"]

    filename = f"{protein_name}_RFD_prediction.pdb"
    (output_dir/filename).write_text(RFD_response, encoding="utf-8")
    # path = os.path.join(output_dir, filename)
    # with open(path, "w", encoding="utf-8") as f:
    #     f.write(RFD_response)

    # visualise RFD results
    html_path = visualize_pdb_py3dmol(RFD_response)
    print("Opened viewer for RF-DIFFUSION:", html_path)

    # ----------------  RUN ProteinMPNN NIM -------------------
    print("Running Protein-MPNN...")
    # use the output of RFD (either access via file, or directly)
    #PMPNN_input = RFD_response
    PMPNN_input = Path(f"{output_dir}/{protein_name}_RFD_prediction.pdb")
    PMPNN_input = filter(lambda line: line.startswith("ATOM"), PMPNN_input.read_text().split("\n"))
    PMPNN_input = "\n".join(list(PMPNN_input)[:400])

    PMPNN_response = PMPNN_NIM(API_KEY, PMPNN_input)
    PMPNN_response = json.loads(PMPNN_response)["mfasta"]

    filename = f"{protein_name}_Protein_MPNN_prediction.fa"
    (output_dir/filename).write_text(PMPNN_response, encoding="utf-8")
    # path = os.path.join(output_dir, filename)
    # with open(path, "w", encoding="utf-8") as f:
    #     f.write(PMPNN_response)

    # --------------------  RUN AF-MULTIMER -------------------
    print("Running ALPHAFOLD-MULTIMER...")

    # ADD read PMPNN_response from text file to get fasta_sequences

    # obtain predicted sequences, note: need to skip first 2 lines
    fasta_sequences = [x.strip() for x in PMPNN_response.split("\n") if '>' not in x][2:]
    # create a list of pairs of [possible_binder_seq, target_seq]
    binder_target_pairs = [[binder, example_input] for binder in fasta_sequences]
    print(f"Generated {len(fasta_sequences)} FASTA sequences and {len(binder_target_pairs)} binder-target pairs.")

    multimer_response_codes = [0 for i in binder_target_pairs]
    multimer_results = [None for i in binder_target_pairs]

    n_processed = 0
    pairs_to_process = 1

    for binder_target_pair in binder_target_pairs:
        if n_processed >= pairs_to_process:
            break
        print(f"MULTIMER: Processing pair number {n_processed+1} of {pairs_to_process}")
        multimer_response = MULTIMER_NIM(API_KEY, binder_target_pair)

        # MULTIMER RETURNS A ZIP, save the zip file
        zip_file = f"{protein_name}_MULTIMER_{n_processed+1}.zip"
        zip_path = os.path.join(output_dir, zip_file)
        with open(zip_path, "wb") as f:
            f.write(multimer_response.content)  # raw bytes
        print("Saved ZIP to:", zip_path)

        # unzip the zip file, then decode contents to get pdb text
        # NEED TO TEST WHEN USING MULTIPLE MODELS
        with zipfile.ZipFile(zip_path, 'r') as zf:
            names = zf.namelist()
            for name in names:
                raw = zf.read(name)
                text = raw.decode("utf-8", "ignore")
        
        obj = json.loads(text)
        pdb_text = "".join(obj)

        # # this method directly uses zip file in memory, but will not save zip file
        # with zipfile.ZipFile(io.BytesIO(multimer_response.content), "r") as zf:
        #     response_name = next(n for n in zf.namelist() if n.lower().endswith(".response"))
        #     raw = zf.read(response_name)                        # bytes
        #     obj = json.loads(raw.decode("utf-8", "ignore"))     # typically list[str]
        #     pdb_text = "".join(obj) if isinstance(obj, list) else str(obj)

        multimer_response_codes[n_processed] = multimer_response.status_code
        multimer_results[n_processed] = pdb_text

        print(f"Finished binder-target pair number {n_processed+1} of {pairs_to_process}")
        n_processed += 1
    
    for i, res in enumerate(multimer_results):
        if not res:
            continue

        filename = f"{protein_name}_pdb_{i+1}_MULTIMER.pdb"
        (output_dir/filename).write_text(res, encoding="utf-8")
        
        # path = os.path.join(output_dir, filename)
        # with open(path, "w", encoding="utf-8") as f:
        #     f.write(res)

        html_path = visualize_pdb_py3dmol(res)
        print("Opened viewer for MULTIMER:", html_path)

    # --------------------  CALCULATE pLDDT SCORE -------------------
    print("calculating pLDDT scores")
    plddts = []
    for i in range(0, len(multimer_results)):
        if multimer_results[i]:
            plddts.append(calculate_average_pLDDT(multimer_results[i]))
    
    print(f"FINAL pLDDT SCORE: {plddts}")

    with open (output_dir/ "pLDDT_scores.txt", "w", encoding="utf-8") as f:
        for i, score in enumerate(plddts):
            f.write(f"{protein_name}_pdb_{i+1}_MULTIMER.pdb\t{score:.6f}\n")





        

        












if __name__ == "__main__":
    main()