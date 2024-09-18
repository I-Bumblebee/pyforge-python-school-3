from os import getenv

from fastapi import FastAPI

from routes.molecule_routes import router as molecule_routes

app = FastAPI()

app.include_router(molecule_routes)


@app.get("/")
def get_server():
    return {"server_id": getenv("SERVER_ID", "1")}
