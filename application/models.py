from django.db import models
from datetime import datetime
from django.contrib.auth.models import User, auth
# Create your models here.


class txtFile(models.Model):
    name = models.FileField(upload_to = "allTxtFiles")
    timestamp = models.DateTimeField(default = datetime.now())