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

table: sa.Table = sa.table(
    "molecules",
    sa.column("id"),
    sa.column("identifier"),
    sa.column("smiles"),
)


def upgrade():
    schema_upgrades()
    if context.get_x_argument(as_dictionary=True).get('data', None):
        data_upgrades()


def downgrade():
    if context.get_x_argument(as_dictionary=True).get('data', None):
        data_downgrades()
    schema_downgrades()


def schema_upgrades():
    op.create_table(
        'molecules',
        sa.Column('id', sa.Integer(), primary_key=True),
        sa.Column('identifier', sa.String(), unique=True),
        sa.Column('smiles', sa.String()),
    )


def schema_downgrades():
    op.drop_table('molecules')


def data_upgrades():
    op.bulk_insert(
        table, [
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
    op.execute(table.delete())
