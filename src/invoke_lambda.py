import json

import boto3

from configs.settings import settings

session = boto3.Session(
    aws_access_key_id=settings.aws_access_key_id,
    aws_secret_access_key=settings.aws_secret_access_key,
    region_name=settings.region_name,
)

lambda_client = session.client('lambda')

event = {"name": "Student"}

response = lambda_client.invoke(
    FunctionName='HelloStudentFunction',
    InvocationType='RequestResponse',
    Payload=json.dumps(event)
)

response_payload = response['Payload'].read()
response_data = json.loads(response_payload)

print(response_data)
