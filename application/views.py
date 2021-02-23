from django.shortcuts import render, redirect
from .models import *
from .RanGen import *
from .InpCode import *
from django.contrib import messages
from django.http import HttpResponse
from django.http import JsonResponse
from django.conf import settings
import json
import os


def tables(request):
    return render(request, "table.html")

# main page function
def main(request):
    return render(request,'main.html')

def index(request):
    context = {}
    # if not request.user.is_authenticated:
    #     return redirect("login")
    if request.method == "POST":
        if "form-steps" in request.POST:
            form_step = request.POST.get("form-steps")
            if form_step and form_step != "":   
                if form_step == "specifying-parameters":
                    n = int(request.POST.get("n"))
                    k = int(request.POST.get("k"))
                    d = int(request.POST.get("d"))
                    q = int(request.POST.get("q"))
                    context['response'] = alltogether(n, k, d, q)
                    context["specifyingParameters"] = "checked" if (form_step == "specifying-parameters") else ""
                
                elif form_step == "input-code":
                    myTxtFile = request.FILES.get("myTxtFile")

                    new_file = txtFile(name = myTxtFile)
                    new_file.save()

                    ip = int(request.POST.get("ip"))
                    iq = int(request.POST.get("iq"))

                    # this is the selected radio option from 
                    # 1. concial
                    # 2. invariant
                    # 3. model p code
                    form_radio = request.POST.get("form-radio")
                    print(form_radio)
                    
                    context["inputCode"] = "checked" if (form_step == "input-code") else ""

                    project_path = settings.BASE_DIR
                    file_path = os.path.join(project_path, "media")
                    file_path = os.path.join(file_path, str(new_file.name))

                    with open(file_path) as myFile:
                        file_content = myFile.read().splitlines()
                        file_content = json.loads(" ".join(file_content))
                        context['file_content'] = file_content
                        if form_radio == "canonical":
                            answers = alltogetherCan(file_content, iq, ip)
                        elif form_radio == "Invariant":
                            answers = alltogetherInv(file_content, iq, ip)
                        else:
                            answers = alltogetherModP(file_content, iq, ip)
                        context['response'] = answers[0]
                        context['distance'] = answers[1]
                        context['preserved'] = answers[2]
                        print(file_content)

                    myFile.close()


    return render(request, 'main.html', context)

# function for signup

def signup(request):
    if request.method == "POST":
        name = request.POST['name']
        l_name = request.POST['l_name']
        email = request.POST['email']
        pass1 = request.POST['pass1']
        pass2 = request.POST['pass2']
        context = {
            "name":name,
            "l_name":l_name,
            "email":email,
            "pass1":pass1,
            "pass2":pass2,
        }
        if pass1==pass2:
            if User.objects.filter(username=email).exists():
                print("Email already taken")
                messages.info(request, "Entered email already in use!")
                context['border'] = "email" 
                return render(request, "signup.html", context)

            user = User.objects.create_user(username=email, first_name=name, password=pass1, last_name=l_name)
            user.save()
            
            return redirect("login")
        else:
            messages.info(request, "Your pasword doesn't match!")
            context['border'] = "password"
            return render(request, "signup.html", context)


    
    return render(request, "signup.html")


# function for login

def login(request):

    if request.method == "POST":
        email = request.POST['email']
        password = request.POST['password']
        context = {
            'email': email,
            'password': password
        }
        user = auth.authenticate(username=email, password=password)
        if user is not None:
            auth.login(request, user)
            return redirect("index")
        else:
            messages.info(request, "Incorrect login details!")
            return render(request, "login.html", context)
            # return redirect("login")
    else:
        return render(request, "login.html")


# function for logout

def logout(request):
    auth.logout(request)
    return redirect("index")

