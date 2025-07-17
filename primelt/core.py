import pandas as pd
import os

# Import your core calculation modules (must be in the same folder)
from . import PRIMELT3P_core
from . import PRIMARSMELT_core

def run_primelt3(input_file, mgo_s=None, feo_s=None, output_dir=None):
    # Set defaults if not provided
    mgo_s = mgo_s if mgo_s is not None else 38.12
    feo_s = feo_s if feo_s is not None else 8.02
    PRIMELT3P_core.config_manager.update_settings(mgo_s=mgo_s, feo_s=feo_s)
    samples = PRIMELT3P_core.read_data(input_file)
    #if not PRIMELT3P_core.test_data(samples):
    #    raise ValueError("Input data format is invalid.")
    solution = PRIMELT3P_core.work(samples)
    magma_afm, magma_batch = solution[1], solution[0]

    # Save results
    base = os.path.splitext(os.path.basename(input_file))[0]
    output_dir = output_dir or os.path.dirname(input_file)
    magma_afm.to_csv(os.path.join(output_dir, f"{base}_afm.csv"), sep=',')
    magma_batch.to_csv(os.path.join(output_dir, f"{base}_batch.csv"), sep=',')

    # Optionally, save Excel
    try:
        with pd.ExcelWriter(os.path.join(output_dir, f"{base}_all.xlsx")) as writer:
            samples.to_excel(writer, sheet_name='original_comp', index=False)
            magma_batch.to_excel(writer, sheet_name='batch')
            magma_afm.to_excel(writer, sheet_name='afm')
    except Exception as e:
        print("Excel file could not be generated:", e)

    print("Results saved to", output_dir)
    return magma_afm, magma_batch

def run_primarsmelt(input_file, mgo_s=None, feo_s=None, output_dir=None):
    # Set defaults if not provided
    mgo_s = mgo_s if mgo_s is not None else 30.2
    feo_s = feo_s if feo_s is not None else 17.9
    PRIMARSMELT_core.config_manager.update_settings(mgo_s=mgo_s, feo_s=feo_s)
    samples = PRIMARSMELT_core.read_data(input_file)
    if not PRIMARSMELT_core.test_data(samples):
        raise ValueError("Input data format is invalid.")
    solution = PRIMARSMELT_core.final_results(samples)
    magma_afm, magma_batch = solution[1], solution[0]

    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    # Save results
    base = os.path.splitext(os.path.basename(input_file))[0]
    output_dir = output_dir or os.path.dirname(input_file)
    magma_afm.to_csv(os.path.join(output_dir, f"{base}_afm.csv"), sep=',')
    magma_batch.to_csv(os.path.join(output_dir, f"{base}_batch.csv"), sep=',')

    # Optionally, save Excel
    try:
        with pd.ExcelWriter(os.path.join(output_dir, f"{base}_all.xlsx")) as writer:
            samples.to_excel(writer, sheet_name='original_comp', index=False)
            magma_batch.to_excel(writer, sheet_name='batch')
            magma_afm.to_excel(writer, sheet_name='afm')
    except Exception as e:
        print("Excel file could not be generated:", e)



    print("Results saved to", output_dir)
    return magma_afm, magma_batch