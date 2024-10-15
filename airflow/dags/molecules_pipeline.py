import io
from datetime import datetime
from typing import Dict, List

import boto3
import pandas as pd
from airflow.decorators import dag, task
from airflow.utils.dates import days_ago
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors

from configs.database import get_db_session
from configs.settings import settings
from src.dao.molecule_dao import MoleculeDAO
from src.models.molecule import Molecule
from src.utils import is_lipinski_pass


@dag(
    dag_id='molecules_pipeline',
    schedule_interval='@daily',
    start_date=days_ago(1),
    catchup=False,
)
def molecules_pipeline():
    @task
    async def extract_data() -> List[Molecule]:
        async with get_db_session() as session:
            molecules_dao = MoleculeDAO(session)
            return [
                molecule async for molecule in molecules_dao.list_molecules()
            ]

    @task
    def transform_data(molecules: List[Molecule]) -> List[Dict]:
        transformed_data = []
        for molecule in molecules:
            mol = Chem.MolFromSmiles(molecule.smiles)
            transformed_data.append(
                {
                    'id': molecule.id,
                    'identifier': molecule.identifier,
                    'smiles': molecule.smiles,
                    'molecular_weight': Descriptors.ExactMolWt(mol),
                    'logP': Crippen.MolLogP(mol),
                    'TPSA': Descriptors.TPSA(mol),
                    'H_donors': Descriptors.NumHDonors(mol),
                    'H_acceptors': Descriptors.NumHAcceptors(mol),
                    'lipinski_pass': is_lipinski_pass(mol)
                }
            )

        return transformed_data

    @task
    def save_transformed_data(transformed_data: List[Dict]):
        df = pd.DataFrame(transformed_data)

        with io.BytesIO() as output:
            with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
                df.to_excel(writer, index=False, sheet_name='Molecules')
            data = output.getvalue()

        s3 = boto3.resource(
            "s3",
            aws_access_key_id=settings.aws_access_key_id,
            aws_secret_access_key=settings.aws_secret_access_key,
            region_name=settings.region_name,
        )

        current_date = datetime.now().strftime('%d_%m_%Y')

        s3.Bucket('hw-bucket-lt').put_object(
            Key=f"transformed_molecules_{current_date}.xlsx",
            Body=data,
        )

    molecules = extract_data()
    transformed_data = transform_data(molecules)
    save_transformed_data(transformed_data)


molecules_dag = molecules_pipeline()
