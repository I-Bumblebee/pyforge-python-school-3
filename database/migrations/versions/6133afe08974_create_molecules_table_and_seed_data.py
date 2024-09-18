"""create_molecules_table_and_seed_data

Revision ID: 6133afe08974
Revises: None
Create Date: 2024-08-26 20:36:32.806081

"""

# revision identifiers, used by Alembic.
revision = '6133afe08974'
down_revision = None

from alembic import op
import sqlalchemy as sa

from alembic import context

molecules_table: sa.Table = None

def upgrade():
    schema_upgrades()
    if context.get_x_argument(as_dictionary=True).get('data', None):
        data_upgrades()

def downgrade():
    if context.get_x_argument(as_dictionary=True).get('data', None):
        data_downgrades()
    schema_downgrades()

def schema_upgrades():
    global molecules_table
    molecules_table = op.create_table('molecules',
        sa.Column('id', sa.Integer()),
        sa.Column('identifier', sa.String(), nullable=True),
        sa.Column('smiles', sa.Text()),
        sa.PrimaryKeyConstraint('id')
    )
    op.create_index(op.f('ix_molecules_id'), 'molecules', ['id'])
    op.create_index(op.f('ix_molecules_identifier'), 'molecules', ['identifier'], unique=True)

def schema_downgrades():
    op.drop_table('molecules')

def data_upgrades():
    global molecules_table
    op.bulk_insert(molecules_table,
        [
            {
                'identifier': 'mol1',
                'smiles': 'C1=CC=CC=C1'
            },
            {
                'identifier': 'mol2',
                'smiles': 'C2H6'
            },
        ]               
    ) 

def data_downgrades():
    global molecules_table
    op.execute(molecules_table.delete())
