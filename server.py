'''

What is this?
server.py is a very dummy static server code only for testing.
It basically just lets you view files in the same folder using localhost:80.
Feel free to change port numbers and stuff.
Make sure we don't send server.py when we update the website :sob:

'''

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles

app = FastAPI()

origins = ["*"]

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

app.mount('/', StaticFiles(directory='./'))

if __name__ == '__main__':
    import uvicorn
    uvicorn.run(app, host='0.0.0.0', port=80)
